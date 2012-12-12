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
  const SingleParticleOrbitals& M_init,
  const Jastrow& v_init,
  cl_uint N_init,
  cl_uint update_hop_maxdist_init,
  const vector<cl_fptype>& t_init,
  cl_fptype U_init,
  cl_fptype W_deviation_target_init,
  cl_uint updates_until_W_recalc_init,
  cl_fptype T_deviation_target_init,
  cl_uint updates_until_T_recalc_init )
  : rng( rng_init ),
    lat( lat_init ), M( M_init ), v( v_init ),
    update_hop_maxdist( update_hop_maxdist_init ),
    t( t_init ), U( U_init ),
    econf( ElectronConfiguration( lat, N_init, &rng ) ),
    Wbu_1(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, N_init )
    ),
    Wbu_2(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, N_init )
    ),
    Wbu_active(   &Wbu_1 ),
    Wbu_inactive( &Wbu_2 ),
    Wd_1(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    Wd_2(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    Wd_active(   M_init.ssym ? &Wd_1 : nullptr ),
    Wd_inactive( M_init.ssym ? &Wd_2 : nullptr ),
    T( Eigen::VectorXfp( lat->L ) ),
    completed_mcsteps( 0 ),
    updates_until_W_recalc( updates_until_W_recalc_init ),
    updates_until_T_recalc( updates_until_T_recalc_init ),
    updates_since_W_recalc( 0 ), updates_since_T_recalc( 0 ),
    W_devstat( FPDevStat( W_deviation_target_init ) ),
    T_devstat( FPDevStat( T_deviation_target_init ) )
{

  do { // loop to get a large overlap

    bool invertible;

    do { // loop to get any overlap at all

      econf.distribute_random();

      if ( M.ssym == true ) {

        Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Du().transpose() );
        invertible = ( lu_decomp.rank() == lu_decomp.rows() );
        if ( invertible ) {
          Wbu_active->noalias()
            = lu_decomp.solve( M.orbitals.transpose() ).transpose();
        } else {
          continue;
        }

        lu_decomp.compute( calc_Dd().transpose() );
        invertible &= ( lu_decomp.rank() == lu_decomp.rows() );
        if ( invertible ) {
          Wd_active ->noalias()
            = lu_decomp.solve( M.orbitals.transpose() ).transpose();
        }

      } else {

        Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Db() );
        invertible = ( lu_decomp.rank() == lu_decomp.rows() );
        if ( invertible ) {
          Wbu_active->noalias()
            = lu_decomp.solve( M.orbitals.transpose() ).transpose();
        }

      }

#if VERBOSE >= 1
      cout << "HubbardModelVMC::HubbardModelVMC() : matrix D is not invertible!"
           << endl;
#endif

    } while ( !invertible );
    // initialize the electrons so that D is invertible
    // (there must be a non-zero overlap between the slater det and |x>)

    // calculate the vector T from scratch
    T = calc_new_T();

    // repeat everything if the initial state has a very low overlap
    // (it's bad for floating point precision before the first recalc)
  } while (

    Wbu_active->array().square().sum()
    / static_cast<cl_fptype>( Wbu_active->size() ) > 10.f ||

    (
      M.ssym == true ?

      Wd_active->array().square().sum()
      / static_cast<cl_fptype>( Wd_active->size() ) > 10.f :

      false
    ) ||

    T.array().square().sum() / static_cast<cl_fptype>( T.size() ) > 10.f

  );

#if VERBOSE >= 1
  cout << "HubbardModelVMC::HubbardModelVMC() : calculated initial matrix W ="
       << endl << *Wbu_active << endl;
  if ( M.ssym == true ) {
    cout << "----->" << endl << *Wd_active << endl;
  }
#endif
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
  for ( cl_uint s = 0; s < econf.N(); ++s ) {
#if VERBOSE >= 1
    cout << "HubbardModelVMC::mcs() : Monte Carlo step = " << completed_mcsteps
         << ", Metropolis step = " << s << endl;
#endif
    metstep();
  }
  ++completed_mcsteps;
}



void HubbardModelVMC::equilibrate( cl_uint N_mcs_equil )
{
  for ( cl_uint n = 0; n < N_mcs_equil; ++n ) {
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

    const cl_fptype R_j = T( lat->get_spinup_site( phop.l ) )
                       / T( lat->get_spinup_site( phop.k_pos ) )
                       * v.exp( 0, 0 ) / v.exp( phop.l, phop.k_pos );

    const cl_fptype R_s = ( M.ssym == true && phop.k >= econf.N() / 2 ) ?
                       ( *Wd_active  )( phop.l - lat->L, phop.k - econf.N() / 2 ) :
                       ( *Wbu_active )( phop.l, phop.k );
    const cl_fptype accept_prob = R_j * R_j * R_s * R_s;

#if VERBOSE >= 1
    cout << "HubbardModelVMC::metstep() : hop possible -> "
         << "R_j = " << R_j
         << ", sdwf_ratio = " << R_s
         << ", accept_prob = " << accept_prob << endl;
#endif

    if ( accept_prob >= 1.f ||
         uniform_real_distribution<cl_fptype>( 0.f, 1.f )( rng ) < accept_prob ) {

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
    if ( M.ssym == true && hop.k >= econf.N() / 2 ) {
      calc_qupdated_Wd( hop );
    } else {
      calc_qupdated_Wbu( hop );
    }

    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new_W();

    cl_fptype dev = calc_deviation( *Wbu_inactive, *Wbu_active );
    if ( M.ssym == true ) {
      dev += calc_deviation( *Wd_inactive, *Wd_active );
    }
    W_devstat.add( dev );

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : recalculated W "
         << "with deviation = " << dev << endl;

    if ( dev > W_devstat.target ) {
      cout << "HubbardModelVMC::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "HubbardModelVMC::perform_W_update() : approximate W =" << endl
           << *Wbu_inactive << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *Wd_inactive << endl;
      }
      cout << "HubbardModelVMC::perform_W_update() : exact W =" << endl
           << *Wbu_active << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *Wd_active << endl;
      }
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
    if ( M.ssym == true && hop.k >= econf.N() / 2 ) {
      calc_qupdated_Wd( hop );
    } else {
      calc_qupdated_Wbu( hop );
    }

#ifndef NDEBUG
    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new_W();

    // swap the buffers (since we want the updated buffer to be the active one)
    swap( Wbu_inactive, Wbu_active );
    if ( M.ssym == true ) {
      swap( Wd_inactive, Wd_active );
    }

    // updated W should now be in the active buffer
    // debug check recalc W should be in the inactive buffer

    cl_fptype dev = calc_deviation( *Wbu_inactive, *Wbu_active );
    if ( M.ssym == true ) {
      dev += calc_deviation( *Wd_inactive, *Wd_active );
    }

# if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > W_devstat.target ) {
      cout << "HubbardModelVMC::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "HubbardModelVMC::perform_W_update() : quickly updated W =" << endl
           << *Wbu_active << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *Wd_active << endl;
      }
      cout << "HubbardModelVMC::perform_W_update() : exact W =" << endl
           << *Wbu_inactive << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *Wd_inactive << endl;
      }
    }
# endif
#endif

    assert( dev < W_devstat.target );
  }
}



Eigen::MatrixXfp HubbardModelVMC::calc_Db() const
{
  assert( M.ssym == false );

  Eigen::MatrixXfp Db( econf.N(), econf.N() );
  for ( cl_uint eid = 0; eid < econf.N(); ++eid ) {
    Db.row( eid ) = M.orbitals.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_Db() : Db = " << endl << Db << endl;
#endif

  return Db;
}



Eigen::MatrixXfp HubbardModelVMC::calc_Du() const
{
  assert( M.ssym == true );

  Eigen::MatrixXfp Du( econf.N() / 2, econf.N() / 2 );
  for ( cl_uint eid = 0; eid < econf.N() / 2; ++eid ) {
    Du.row( eid ) = M.orbitals.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_Du() : Du = " << endl << Du << endl;
#endif

  return Du;
}



Eigen::MatrixXfp HubbardModelVMC::calc_Dd() const
{
  assert( M.ssym == true );

  Eigen::MatrixXfp Dd( econf.N() / 2, econf.N() / 2 );
  for ( cl_uint eid = econf.N() / 2; eid < econf.N(); ++eid ) {
    Dd.row( eid - econf.N() / 2 )
      = M.orbitals.row( lat->get_spinup_site( econf.get_electron_pos( eid ) ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_Dd() : Dd = " << endl << Dd << endl;
#endif

  return Dd;
}


void HubbardModelVMC::calc_new_W()
{
  if ( M.ssym == true ) {

    Wbu_inactive->noalias()
      = calc_Du().transpose().partialPivLu().solve( M.orbitals.transpose() ).transpose();
    Wd_inactive->noalias()
      = calc_Dd().transpose().partialPivLu().solve( M.orbitals.transpose() ).transpose();

    swap( Wbu_inactive, Wbu_active );
    swap( Wd_inactive, Wd_active );

  } else {

    Wbu_inactive->noalias()
      = calc_Db().transpose().partialPivLu().solve( M.orbitals.transpose() ).transpose();

    swap( Wbu_inactive, Wbu_active );

  }
}



void HubbardModelVMC::calc_qupdated_Wbu( const ElectronHop& hop )
{
  *Wbu_inactive = *Wbu_active;

  Wbu_inactive->noalias() -=
    ( Wbu_active->col( hop.k ) / ( *Wbu_active )( hop.l, hop.k ) )
    * ( Wbu_active->row( hop.l ) - Wbu_active->row( hop.k_pos ) );

  swap( Wbu_inactive, Wbu_active );
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

    cl_fptype dev = calc_deviation( T_approx, T );
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
    cl_fptype dev = calc_deviation( T, T_chk );

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

  for ( cl_uint i = 0; i < lat->L; ++i ) {
    cl_fptype sum = 0.f;
    for ( cl_uint j = 0; j < lat->L; ++j ) {
      sum += v( i, j ) * static_cast<cl_fptype>(
               ( econf.get_site_occ( j ) + econf.get_site_occ( j + lat->L ) ) );
    }
    T_new( i ) = exp( sum );
  }

  return T_new;
}



Eigen::VectorXfp HubbardModelVMC::calc_qupdated_T( const ElectronHop& hop ) const
{
  Eigen::VectorXfp T_prime( lat->L );

  for ( cl_uint i = 0; i < lat->L; ++i ) {
    T_prime( i ) = T( i ) * v.exp( i, lat->get_spinup_site( hop.l ) )
                   / v.exp( i, lat->get_spinup_site( hop.k_pos ) );
  }

  return T_prime;
}



cl_fptype HubbardModelVMC::E_l()
{
  // calculate expectation value of the T part of H
  cl_fptype E_l_kin = 0.f;

  for ( cl_uint k = 0; k < econf.N(); ++k ) {

    const cl_uint k_pos = econf.get_electron_pos( k );
    assert( econf.get_site_occ( k_pos ) == ELECTRON_OCCUPATION_FULL );

    for ( cl_uint X = 1; X <= t.size(); ++X ) {
      if ( t[X - 1] == 0.f ) {
        continue;
      }

      cl_fptype sum_Xnn = 0.f;
      lat->get_Xnn( k_pos, X, &k_pos_Xnn );
      for ( auto l_it = k_pos_Xnn.begin(); l_it != k_pos_Xnn.end(); ++l_it ) {
        if ( econf.get_site_occ( *l_it ) == ELECTRON_OCCUPATION_EMPTY ) {
          const cl_fptype R_j =   T( lat->get_spinup_site( *l_it ) )
                               / T( lat->get_spinup_site( k_pos ) )
                               * v.exp_onsite() / v.exp( *l_it, k_pos );
          if ( M.ssym == true && k >= econf.N() / 2 ) {
            sum_Xnn += R_j * ( *Wd_active )( *l_it - lat->L, k - econf.N() / 2 );
          } else {
            sum_Xnn += R_j * ( *Wbu_active )( *l_it, k );
          }
        }
      }
      E_l_kin -= t[X - 1] * sum_Xnn;

    }
  }

  const cl_fptype E_l_result =
    ( E_l_kin + U * econf.get_num_dblocc() ) /
    static_cast<cl_fptype>( lat->L );

#if VERBOSE >= 1
  cout << "HubbardModelVMC::E_l() = " << E_l_result << endl;
#endif

  return E_l_result;
}



cl_ulong HubbardModelVMC::mctime() const
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
