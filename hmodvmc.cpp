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

#if VERBOSE >= 1
# include <iostream>
#endif

#include <cmath>
#include <algorithm>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/LU>

#ifdef USE_CBLAS
extern "C" {
# include <cblas.h>
}
#endif

using namespace std;



HubbardModelVMC::HubbardModelVMC(
  const shared_ptr<mt19937>& rng_init,
  const shared_ptr<Lattice>& lat_init,
  const SingleParticleOrbitals& detwf_init,
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
    lat( lat_init ), detwf( detwf_init ), v( v_init ),
    update_hop_maxdist( update_hop_maxdist_init ),
    t( t_init ), U( U_init ),
    econf( ElectronConfiguration( lat, N_init, rng ) ),
    Wbu_1(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, N_init )
    ),
    Wbu_2(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, N_init )
    ),
    Wbu_active(   &Wbu_1 ),
    Wbu_inactive( &Wbu_2 ),
    Wd_1(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    Wd_2(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    Wd_active(   detwf_init.ssym ? &Wd_1 : nullptr ),
    Wd_inactive( detwf_init.ssym ? &Wd_2 : nullptr ),
    T( Eigen::VectorXfp( lat->L ) ),
#ifdef USE_CBLAS
    tempWcol(
      detwf_init.ssym ?
      Eigen::VectorXfp( lat->L ) :
      Eigen::VectorXfp( 2 * lat->L )
    ),
    tempWrow(
      detwf_init.ssym ?
      Eigen::VectorXfp( N_init / 2 ) :
      Eigen::VectorXfp( N_init )
    ),
#endif
    updates_until_W_recalc( updates_until_W_recalc_init ),
    updates_until_T_recalc( updates_until_T_recalc_init ),
    updates_since_W_recalc( 0 ), updates_since_T_recalc( 0 ),
    W_devstat( FPDevStat( W_deviation_target_init ) ),
    T_devstat( FPDevStat( T_deviation_target_init ) )
{

  bool enough_overlap;

  do {

    econf.distribute_random();
    enough_overlap = true; // assume until we are proven wrong

#if VERBOSE >=1
    cout << "HubbardModelVMC::HubbardModelVMC() : checking newly generated "
         << "state for enough overlap" << endl;
#endif

    fptype Wbu_avg, Wd_avg;

    if ( detwf.ssym == true ) {

      // check spin up part
      Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Du().transpose() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : spin up part has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      Wbu_active->noalias()
        = lu_decomp.solve( detwf.M.transpose() ).transpose();
      Wbu_avg
        = Wbu_active->squaredNorm() / static_cast<fptype>( Wbu_active->size() );
      enough_overlap &= Wbu_avg < 100.f;
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : spin up part has too "
             << "little overlap with the determinantal wavefunction, "
             << "inverse overlap measure is: " << Wbu_avg << endl;
#endif
        continue;
      }

      // check spin down part
      lu_decomp.compute( calc_Dd().transpose() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : spin down part has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      Wd_active->noalias()
        = lu_decomp.solve( detwf.M.transpose() ).transpose();
      Wd_avg
        = Wd_active->squaredNorm()  / static_cast<fptype>( Wd_active->size() );
      enough_overlap &= Wd_avg < 100.f;
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : spin down part has too "
             << "little overlap with the determinantal wavefunction, "
             << "inverse overlap measure is: " << Wd_avg << endl;
#endif
        continue;
      }

    } else {

      // check whole determinantal part
      Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Db().transpose() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : state has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      Wbu_active->noalias()
        = lu_decomp.solve( detwf.M.transpose() ).transpose();
      Wbu_avg
        = Wbu_active->squaredNorm() / static_cast<fptype>( Wbu_active->size() );
      enough_overlap &= Wbu_avg < 50.f;
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : state has too "
             << "little overlap with the determinantal wavefunction, "
             << "inverse overlap measure is: " << Wbu_avg << endl;
#endif
        continue;
      }

    }

    // check Jastrow part
    T = calc_new_T();
    fptype T_avg = T.squaredNorm() / static_cast<fptype>( T.size() );
    enough_overlap &= T_avg < 100.f;
    if ( !enough_overlap ) {
#if VERBOSE >= 1
      cout << "HubbardModelVMC::HubbardModelVMC() : Jastrow ratios "
           << "are to small, inverse measure is: " << T_avg << endl;
#endif
      continue;
    }

#if VERBOSE >= 1
    cout << "HubbardModelVMC::HubbardModelVMC() : state has sufficient "
         << "overlap! inverse overlap measures are: "
         << Wbu_avg << " " << Wd_avg << " " << T_avg
         << " -> initial state selection completed!" << endl;
#endif

  } while ( !enough_overlap );

#if VERBOSE >= 1
  cout << "HubbardModelVMC::HubbardModelVMC() : calculated initial matrix W ="
       << endl << *Wbu_active << endl;
  if ( detwf.ssym == true ) {
    cout << "----->" << endl << *Wd_active << endl;
  }
#endif

}



void HubbardModelVMC::mcs()
{
#if VERBOSE >= 1
  cout << "HubbardModelVMC::mcs() : starting new Monte Carlo step!" << endl;
#endif

  // perform a number of metropolis steps equal to the number of electrons
  for ( unsigned int s = 0; s < lat->L; ++s ) {
#if VERBOSE >= 1
    cout << "HubbardModelVMC::mcs() : Metropolis step = " << s << endl;
#endif
    metstep();
  }
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

    const fptype R_j = T( lat->get_spinup_site( phop.l ) )
                       / T( lat->get_spinup_site( phop.k_pos ) )
                       * v.exp_onsite() / v.exp( phop.l, phop.k_pos );

    const fptype R_s = ( detwf.ssym == true && phop.k >= econf.N() / 2 ) ?
                       ( *Wd_active  )( phop.l - lat->L, phop.k - econf.N() / 2 ) :
                       ( *Wbu_active )( phop.l, phop.k );
    const fptype accept_prob = R_j * R_j * R_s * R_s;

#if VERBOSE >= 1
    cout << "HubbardModelVMC::metstep() : hop possible -> "
         << "R_j = " << R_j
         << ", sdwf_ratio = " << R_s
         << ", accept_prob = " << accept_prob << endl;
#endif

    if ( accept_prob >= 1.f ||
         uniform_real_distribution<fptype>( 0.f, 1.f )( *rng ) < accept_prob ) {

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
    calc_qupdated_W( hop );

    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new_W();

    fptype dev = calc_deviation( *Wbu_inactive, *Wbu_active );
    if ( detwf.ssym == true ) {
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
      if ( detwf.ssym == true ) {
        cout << "----->" << endl << *Wd_inactive << endl;
      }
      cout << "HubbardModelVMC::perform_W_update() : exact W =" << endl
           << *Wbu_active << endl;
      if ( detwf.ssym == true ) {
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
    calc_qupdated_W( hop );

#ifndef NDEBUG
    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new_W();

    // swap the buffers (since we want the updated buffer to be the active one)
    swap( Wbu_inactive, Wbu_active );
    if ( detwf.ssym == true ) {
      swap( Wd_inactive, Wd_active );
    }

    // updated W should now be in the active buffer
    // debug check recalc W should be in the inactive buffer

    fptype dev = calc_deviation( *Wbu_inactive, *Wbu_active );
    if ( detwf.ssym == true ) {
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
      if ( detwf.ssym == true ) {
        cout << "----->" << endl << *Wd_active << endl;
      }
      cout << "HubbardModelVMC::perform_W_update() : exact W =" << endl
           << *Wbu_inactive << endl;
      if ( detwf.ssym == true ) {
        cout << "----->" << endl << *Wd_inactive << endl;
      }
    }
# endif
    assert( dev < W_devstat.target );
#endif
  }
}



Eigen::MatrixXfp HubbardModelVMC::calc_Db() const
{
  assert( detwf.ssym == false );

  Eigen::MatrixXfp Db( econf.N(), econf.N() );
  for ( unsigned int eid = 0; eid < econf.N(); ++eid ) {
    Db.row( eid ) = detwf.M.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_Db() : Db = " << endl << Db << endl;
#endif

  return Db;
}



Eigen::MatrixXfp HubbardModelVMC::calc_Du() const
{
  assert( detwf.ssym == true );

  Eigen::MatrixXfp Du( econf.N() / 2, econf.N() / 2 );
  for ( unsigned int eid = 0; eid < econf.N() / 2; ++eid ) {
    Du.row( eid ) = detwf.M.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_Du() : Du = " << endl << Du << endl;
#endif

  return Du;
}



Eigen::MatrixXfp HubbardModelVMC::calc_Dd() const
{
  assert( detwf.ssym == true );

  Eigen::MatrixXfp Dd( econf.N() / 2, econf.N() / 2 );
  for ( unsigned int eid = econf.N() / 2; eid < econf.N(); ++eid ) {
    Dd.row( eid - econf.N() / 2 )
      = detwf.M.row( lat->get_spinup_site( econf.get_electron_pos( eid ) ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_Dd() : Dd = " << endl << Dd << endl;
#endif

  return Dd;
}


void HubbardModelVMC::calc_new_W()
{
  if ( detwf.ssym == true ) {

    Wbu_inactive->noalias()
      = calc_Du().transpose().partialPivLu().solve( detwf.M.transpose() ).transpose();
    Wd_inactive->noalias()
      = calc_Dd().transpose().partialPivLu().solve( detwf.M.transpose() ).transpose();

    swap( Wbu_inactive, Wbu_active );
    swap( Wd_inactive, Wd_active );

  } else {

    Wbu_inactive->noalias()
      = calc_Db().transpose().partialPivLu().solve( detwf.M.transpose() ).transpose();

    swap( Wbu_inactive, Wbu_active );

  }
}



void HubbardModelVMC::calc_qupdated_W( const ElectronHop& hop )
{
  unsigned int k     = hop.k;
  unsigned int l     = hop.l;
  unsigned int k_pos = hop.k_pos;

  Eigen::MatrixXfp*& W = ( detwf.ssym == true && hop.k >= econf.N() / 2 ) ?
                         Wd_active : Wbu_active;

  if ( detwf.ssym == true && hop.k >= econf.N() / 2 ) {
    k     -= econf.N() / 2;
    l     -= lat->L;
    k_pos -= lat->L;
  }

#ifdef USE_CBLAS

  tempWcol = W->col( k );
  tempWrow = W->row( l ) - W->row( k_pos );

#ifdef USE_FP_DBLPREC
  cblas_dger(
#else
  cblas_sger(
#endif
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
    CblasRowMajor,
#else
    CblasColMajor,
#endif
    W->rows(),
    W->cols(),
    - 1.f / ( *W )( l, k ),
    tempWcol.data(),
    1,
    tempWrow.data(),
    1,
    W->data(),
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
    W->cols()
#else
    W->rows()
#endif
  );

#else // #ifndef USE_CBLAS

  Eigen::MatrixXfp*& W_inactive = ( detwf.ssym == true && hop.k >= econf.N() / 2 ) ?
                                  Wd_inactive : Wbu_inactive;

  *W_inactive = *W;

  W_inactive->noalias() -=
    ( W->col( k ) / ( *W )( l, k ) )
    * ( W->row( l ) - W->row( k_pos ) );

  swap( W_inactive, W );

#endif
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
    T_prime( i ) = T( i ) * v.exp( i, lat->get_spinup_site( hop.l ) )
                   / v.exp( i, lat->get_spinup_site( hop.k_pos ) );
  }

  return T_prime;
}



fptype HubbardModelVMC::E_l()
{
  // calculate expectation value of the T part of H
  fptype E_l_kin = 0.f;

  // loop over different elektrons k
  for ( unsigned int k = 0; k < econf.N(); ++k ) {

    const unsigned int k_pos = econf.get_electron_pos( k );
    assert( econf.get_site_occ( k_pos ) == ELECTRON_OCCUPATION_FULL );

    // loop over different neighbor orders X
    for ( unsigned int X = 1; X <= t.size(); ++X ) {
      if ( t[X - 1] == 0.f ) {
        continue;
      }

      fptype sum_Xnn = 0.f;
      lat->get_Xnn( k_pos, X, &k_pos_Xnn );

      // calculate part of R_j that is constant for this X and k
      assert( k_pos_Xnn.size() != 0 );
      const fptype R_j_constXk =
        v.exp_onsite() / v.exp( k_pos_Xnn[0], k_pos )
        / T( lat->get_spinup_site( k_pos ) );
      // (it is possible to do the idxrel reduction only for one of the
      // neighbours as it is guaranteed to be the same for all of them)

      // loop over different neighbours l of order X
      for ( auto l_it = k_pos_Xnn.begin(); l_it != k_pos_Xnn.end(); ++l_it ) {
        if ( econf.get_site_occ( *l_it ) == ELECTRON_OCCUPATION_EMPTY ) {
          const fptype R_j = T( lat->get_spinup_site( *l_it ) ) * R_j_constXk;
          if ( detwf.ssym == true && k >= econf.N() / 2 ) {
            sum_Xnn += R_j * ( *Wd_active )( *l_it - lat->L, k - econf.N() / 2 );
          } else {
            sum_Xnn += R_j * ( *Wbu_active )( *l_it, k );
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



FPDevStat HubbardModelVMC::get_W_devstat() const
{
  return W_devstat;
}
FPDevStat HubbardModelVMC::get_T_devstat() const
{
  return T_devstat;
}
