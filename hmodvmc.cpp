/*
 * Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
 *
 * This file is part of hVMC.
 *
 * hVMC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hVMC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hVMC.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hmodvmc.hpp"
using namespace std;



HubbardModelVMC::HubbardModelVMC(
  mt19937 rng_init,
  Lattice* const lat_init,
  const Eigen::MatrixXfp M_init,
  const Jastrow& v_init,
  unsigned int N_init,
  unsigned int update_hop_maxdist_init,
  const vector<fptype>& t_init,
  fptype U_init,
  fptype W_deviation_target_init,
  unsigned int updates_until_W_recalc_init,
  fptype T_deviation_target_init,
  unsigned int updates_until_T_recalc_init )
  : rng( rng_init ),
    lat( lat_init ), M( M_init ), v( v_init ),
    update_hop_maxdist( update_hop_maxdist_init ),
    t( t_init ), U( U_init ),
    econf( ElectronConfiguration( lat, N_init, &rng ) ),
    Wu_1( Eigen::MatrixXfp( lat->L, N_init / 2 ) ),
    Wu_2( Eigen::MatrixXfp( lat->L, N_init / 2 ) ),
    Wu_active(   &Wu_1 ),
    Wu_inactive( &Wu_2 ),
    Wd_1( Eigen::MatrixXfp( lat->L, N_init / 2 ) ),
    Wd_2( Eigen::MatrixXfp( lat->L, N_init / 2 ) ),
    Wd_active(   &Wd_1 ),
    Wd_inactive( &Wd_2 ),
    T( Eigen::VectorXfp( lat->L ) ),
    completed_mcsteps( 0 ),
    updates_until_W_recalc( updates_until_W_recalc_init ),
    updates_until_T_recalc( updates_until_T_recalc_init ),
    updates_since_W_recalc( 0 ), updates_since_T_recalc( 0 ),
    W_devstat( FPDevStat( W_deviation_target_init ) ),
    T_devstat( FPDevStat( T_deviation_target_init ) )
{
  Eigen::MatrixXfp W( 2 * lat->L, econf.N() );
  do {
    econf.distribute_random();
    Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_D() );

    while ( lu_decomp.isInvertible() == false ) {
      // initialize the electrons so that D is invertible
      // (there must be a non-zero overlap between the slater det and |x>)

#if VERBOSE >= 1
      cout << "HubbardModelVMC::HubbardModelVMC() : matrix D is not invertible!"
           << endl;
#endif
      econf.distribute_random();
      lu_decomp.compute( calc_D() );
    }

    // calculate the W matrix from scratch: W = M * D^-1
    W.noalias() = M * lu_decomp.inverse();

    // calculate the vector T from scratch
    T = calc_new_T();

    // repeat everything if the initial state has a very low overlap
    // (it's bad for floating point precision before the first recalc)
  } while ( W.array().square().sum() / static_cast<fptype>( W.size() ) > 10.f ||
            T.array().square().sum() / static_cast<fptype>( T.size() ) > 10.f );

#if VERBOSE >= 1
  cout << "HubbardModelVMC::HubbardModelVMC() : calculated initial "
       << "W = " << endl << W << endl;
#endif

  // copy from the large matrix W into Wu und Wd
  *Wu_active = W.topLeftCorner(     lat->L, econf.N() / 2 );
  *Wd_active = W.bottomRightCorner( lat->L, econf.N() / 2 );
}



HubbardModelVMC::~HubbardModelVMC()
{
  delete lat;
}



void HubbardModelVMC::mcs()
{
#if VERBOSE >= 1
  cout << "HubbardModelVMC::mcs() : starting new Monte Carlo step!" << endl;
#endif

  // perform a number of metropolis steps equal to the number of electrons
  for ( unsigned int s = 0; s < econf.N(); ++s ) {
#if VERBOSE >= 1
    cout << "HubbardModelVMC::mcs() : Monte Carlo step = " << completed_mcsteps
         << ", Metropolis step = " << s << endl;
#endif
    metstep();
  }
  ++completed_mcsteps;
}



void HubbardModelVMC::equilibrate( unsigned int N_mcs_equil )
{
  for ( unsigned int n = 0; n < N_mcs_equil; ++n ) {
    mcs();
  }
  completed_mcsteps -= N_mcs_equil;
}



bool HubbardModelVMC::metstep()
{
  // let the electron configuration propose a random hop
  const ElectronHop& phop = econf.propose_random_hop( update_hop_maxdist );


  // check if the hop is possible (hopto site must be empty)
  if ( phop.possible == false ) {

    // hop is not possible, rejected!
#if VERBOSE >= 1
    cout << "HubbardModelVMC::metstep() : hop impossible!" << endl;
#endif
    return false;

  } else { // hop possible!

    const fptype R_j = T( lat->get_spinup_site( phop.l ) ) /
                       T( lat->get_spinup_site( phop.k_pos ) ) *
                       exp( v( 0, 0 ) - v( phop.l, phop.k_pos ) );

    const fptype R_s = phop.k < econf.N() / 2 ?
                       ( *Wu_active )( phop.l, phop.k ) :
                       ( *Wd_active )( phop.l - lat->L, phop.k - econf.N() / 2 );
    const fptype accept_prob = R_j * R_j * R_s * R_s;

#if VERBOSE >= 1
    cout << "HubbardModelVMC::metstep() : hop possible -> "
         << "R_j = " << R_j
         << ", sdwf_ratio = " << R_s
         << ", accept_prob = " << accept_prob << endl;
#endif

    if ( accept_prob >= 1.f ||
         uniform_real_distribution<fptype>( 0.f, 1.f )( rng ) < accept_prob ) {

#if VERBOSE >= 1
      cout << "HubbardModelVMC::metstep() : hop accepted!" << endl;
#endif

      econf.do_hop( phop );

      perform_W_update( phop );
      perform_T_update( phop );

      return true;

    } else { // hop possible but rejected!

#if VERBOSE >= 1
      cout << "HubbardModelVMC::metstep() : hop rejected!" << endl;
#endif

      return false;
    }
  }
}



void HubbardModelVMC::perform_W_update( const ElectronHop& hop )
{
  if ( updates_since_W_recalc >= updates_until_W_recalc ) {

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : recalculating W!" << endl;
#endif

    updates_since_W_recalc = 0;

    // puts updated W into the active buffer
    // (only buffer of the hopping spin direction is changed)
    if ( hop.k < econf.N() / 2 ) {
      calc_qupdated_Wu( hop );
    } else {
      calc_qupdated_Wd( hop );
    }

    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new_W();

    fptype dev =   calc_deviation( *Wu_inactive, *Wu_active )
                   + calc_deviation( *Wd_inactive, *Wd_active );
    W_devstat.add( dev );

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : recalculated W "
         << "with deviation = " << dev << endl;

    if ( dev > W_devstat.target ) {
      cout << "HubbardModelVMC::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "HubbardModelVMC::perform_W_update() : approximate W =" << endl
           << *Wu_inactive << endl << endl << *Wd_inactive << endl
           << "HubbardModelVMC::perform_W_update() : exact W =" << endl
           << *Wu_active << endl << endl << *Wd_active << endl;
    }
#endif

    assert( dev < W_devstat.target );

  } else {

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : "
         << "performing a quick update of W!" << endl;
#endif

    ++updates_since_W_recalc;

    // puts updated W into the active buffer
    // (only buffer of the hopping spin direction is changed)
    if ( hop.k < econf.N() / 2 ) {
      calc_qupdated_Wu( hop );
    } else {
      calc_qupdated_Wd( hop );
    }

#ifndef NDEBUG
    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new_W();

    // swap the buffers (since we want the updated buffer to be the active one)
    swap( Wu_inactive, Wu_active );
    swap( Wd_inactive, Wd_active );

    // updated W should now be in the active buffer
    // debug check recalc W should be in the inactive buffer

    fptype dev =   calc_deviation( *Wu_active, *Wu_inactive )
                   + calc_deviation( *Wd_active, *Wd_inactive );

# if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > W_devstat.target ) {
      cout << "HubbardModelVMC::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "HubbardModelVMC::perform_W_update() : quickly updated W =" << endl
           << *Wu_active << endl << endl << *Wd_active << endl
           << "HubbardModelVMC::perform_W_update() : exact W =" << endl
           << *Wu_inactive << endl << endl << *Wd_inactive << endl;
    }
# endif
#endif

    assert( dev < W_devstat.target );
  }
}



Eigen::MatrixXfp HubbardModelVMC::calc_D() const
{
  Eigen::MatrixXfp D( econf.N(), econf.N() );
  for ( unsigned int eid = 0; eid < econf.N(); ++eid ) {
    D.row( eid ) = M.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_D() : D = " << endl << D << endl;
#endif

  return D;
}



void HubbardModelVMC::calc_new_W()
{
  // TODO: LU decomp separately for Wu and Wd (Robert Rueger, 2012-11-26 14:58)
  Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_D() );
  assert( lu_decomp.isInvertible() );
  Eigen::MatrixXfp W = M * lu_decomp.inverse();

  // copy from the large matrix W into the inactive buffers
  *Wu_inactive = W.topLeftCorner(     lat->L, econf.N() / 2 );
  *Wd_inactive = W.bottomRightCorner( lat->L, econf.N() / 2 );

  // swap active and inactive buffers
  swap( Wu_inactive, Wu_active );
  swap( Wd_inactive, Wd_active );
}



void HubbardModelVMC::calc_qupdated_Wu( const ElectronHop& hop )
{
  *Wu_inactive = *Wu_active;

  Wu_inactive->noalias() -=
    ( Wu_active->col( hop.k ) / ( *Wu_active )( hop.l, hop.k ) )
    * ( Wu_active->row( hop.l ) - Wu_active->row( hop.k_pos ) );

  swap( Wu_inactive, Wu_active );
}



void HubbardModelVMC::calc_qupdated_Wd( const ElectronHop& hop )
{
  *Wd_inactive = *Wd_active;

  Wd_inactive->noalias() -=
    ( Wd_active->col( hop.k - econf.N() / 2 )
      / ( *Wd_active )( hop.l - lat->L, hop.k - econf.N() / 2 ) )
    * ( Wd_active->row( hop.l - lat->L )
        - Wd_active->row( hop.k_pos - lat->L ) );

  swap( Wd_inactive, Wd_active );
}



void HubbardModelVMC::perform_T_update( const ElectronHop& hop )
{
  if ( updates_since_T_recalc >= updates_until_T_recalc ) {

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_T_update() : recalculating T!" << endl;
#endif

    updates_since_T_recalc = 0;

    const Eigen::MatrixXfp& T_approx = calc_qupdated_T( hop );
    T = calc_new_T();

    fptype dev = calc_deviation( T_approx, T );
    T_devstat.add( dev );

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_T_update() : recalculated T "
         << "with deviation = " << dev << endl;

    if ( dev > T_devstat.target ) {
      cout << "HubbardModelVMC::perform_T_update() : deviation goal for matrix "
           << "T not met!" << endl
           << "HubbardModelVMC::perform_T_update() : approximate T =" << endl
           << T_approx.transpose() << endl
           << "HubbardModelVMC::perform_T_update() : exact T =" << endl
           << T.transpose() << endl;
    }
#endif

    assert( dev < T_devstat.target );

  } else {

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_T_update() : "
         << "performing a quick update of T!" << endl;
#endif

    ++updates_since_T_recalc;

    T = calc_qupdated_T( hop );

#ifndef NDEBUG
    const Eigen::MatrixXfp& T_chk = calc_new_T();
    fptype dev = calc_deviation( T, T_chk );

# if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_T_update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > T_devstat.target ) {
      cout << "HubbardModelVMC::perform_T_update() : deviation goal for matrix "
           << "T not met!" << endl
           << "HubbardModelVMC::perform_T_update() : quickly updated T =" << endl
           << T.transpose() << endl
           << "HubbardModelVMC::perform_T_update() : exact T =" << endl
           << T_chk.transpose() << endl;
    }
# endif
#endif

    assert( dev < T_devstat.target );
  }
}



Eigen::VectorXfp HubbardModelVMC::calc_new_T() const
{
  Eigen::VectorXfp T_new( lat->L );

  for ( unsigned int i = 0; i < lat->L; ++i ) {
    fptype sum = 0.f;
    for ( unsigned int j = 0; j < lat->L; ++j ) {
      sum += v( i, j ) * static_cast<fptype>(
               ( econf.get_site_occ( j ) + econf.get_site_occ( j + lat->L ) ) );
    }
    T_new( i ) = exp( sum );
  }

  return T_new;
}



Eigen::VectorXfp HubbardModelVMC::calc_qupdated_T( const ElectronHop& hop ) const
{
  Eigen::VectorXfp T_prime( lat->L );

  for ( unsigned int i = 0; i < lat->L; ++i ) {
    T_prime( i ) = T( i ) * exp(   v( i, lat->get_spinup_site( hop.l ) )
                                   - v( i, lat->get_spinup_site( hop.k_pos ) ) );
  }

  return T_prime;
}



fptype HubbardModelVMC::E_l() const
{
  // calculate expectation value of the T part of H
  fptype E_l_kin = 0.f;

  for ( unsigned int k = 0; k < econf.N(); ++k ) {

    const unsigned int k_pos = econf.get_electron_pos( k );
    assert( econf.get_site_occ( k_pos ) == ELECTRON_OCCUPATION_FULL );

    for ( unsigned int X = 1; X <= t.size(); ++X ) {
      if ( t[X - 1] == 0.f ) {
        continue;
      }

      fptype sum_Xnn = 0.f;
      const vector<unsigned int>& k_pos_Xnn = lat->get_Xnn( k_pos, X );
      for ( auto l_it = k_pos_Xnn.begin(); l_it != k_pos_Xnn.end(); ++l_it ) {
        if ( econf.get_site_occ( *l_it ) == ELECTRON_OCCUPATION_EMPTY ) {
          const fptype R_j = T( lat->get_spinup_site( *l_it ) )
                             / T( lat->get_spinup_site( k_pos ) ) *
                             exp( v( 0, 0 ) - v( *l_it, k_pos ) );
          if ( k < econf.N() / 2 ) {
            sum_Xnn += R_j * ( *Wu_active )( *l_it, k );
          } else {
            sum_Xnn += R_j * ( *Wd_active )( *l_it - lat->L, k - econf.N() / 2 );
          }
        }
      }
      E_l_kin -= t[X - 1] * sum_Xnn;

    }
  }

  const fptype E_l_result =
    ( E_l_kin + U * econf.get_num_dblocc() ) /
    static_cast<fptype>( lat->L );

#if VERBOSE >= 1
  cout << "HubbardModelVMC::E_l() = " << E_l_result << endl;
#endif

  return E_l_result;
}



unsigned long int HubbardModelVMC::mctime() const
{
  return completed_mcsteps;
}



FPDevStat HubbardModelVMC::get_W_devstat() const
{
  return W_devstat;
}
FPDevStat HubbardModelVMC::get_T_devstat() const
{
  return T_devstat;
}
