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

#include <eigen3/Eigen/LU>

#ifdef USE_ATLAS
extern "C" {
# include <cblas.h>
# include <clapack.h>
}
#endif

using namespace std;



HubbardModelVMC::HubbardModelVMC(
  mt19937 rng_init,
  Lattice* const lat_init,
  const SingleParticleOrbitals& M_init,
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
    WbuT_1(
      M_init.ssym ?
      Eigen::MatrixXfp( N_init / 2, lat->L ) :
      Eigen::MatrixXfp( N_init, 2 * lat->L )
    ),
    WbuT_2(
      M_init.ssym ?
      Eigen::MatrixXfp( N_init / 2, lat->L ) :
      Eigen::MatrixXfp( N_init, 2 * lat->L )
    ),
    WbuT_active(   &WbuT_1 ),
    WbuT_inactive( &WbuT_2 ),
    WdT_1(
      M_init.ssym ?
      Eigen::MatrixXfp( N_init / 2, lat->L ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    WdT_2(
      M_init.ssym ?
      Eigen::MatrixXfp( N_init / 2, lat->L ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    WdT_active(   M_init.ssym ? &WdT_1 : nullptr ),
    WdT_inactive( M_init.ssym ? &WdT_2 : nullptr ),
    T( Eigen::VectorXfp( lat->L ) ),
#ifdef USE_ATLAS
    tempWTrow(
      M_init.ssym ?
      Eigen::VectorXfp( lat->L ) :
      Eigen::VectorXfp( 2 * lat->L )
    ),
    tempWTcol(
      M_init.ssym ?
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

    if ( M.ssym == true ) {

      // check spin up part
      Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_DuT() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : spin up part has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      *WbuT_active = lu_decomp.solve( M.orbitalsT );
      Wbu_avg
        = WbuT_active->squaredNorm() / static_cast<fptype>( WbuT_active->size() );
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
      lu_decomp.compute( calc_DdT() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : spin down part has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      *WdT_active = lu_decomp.solve( M.orbitalsT );
      Wd_avg
        = WdT_active->squaredNorm() / static_cast<fptype>( WdT_active->size() );
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
      Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_DbT() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC::HubbardModelVMC() : state has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      *WbuT_active = lu_decomp.solve( M.orbitalsT );
      Wbu_avg
        = WbuT_active->squaredNorm() / static_cast<fptype>( WbuT_active->size() );
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
  cout << "HubbardModelVMC::HubbardModelVMC() : calculated initial matrix W^T ="
       << endl << *WbuT_active << endl;
  if ( M.ssym == true ) {
    cout << "----->" << endl << *WdT_active << endl;
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

    const fptype R_s = ( M.ssym == true && phop.k >= econf.N() / 2 ) ?
                       ( *WdT_active  )( phop.k - econf.N() / 2, phop.l - lat->L ) :
                       ( *WbuT_active )( phop.k, phop.l );
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
    calc_qupdated_W( hop );

    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new_W();

    fptype dev = calc_deviation( *WbuT_inactive, *WbuT_active );
    if ( M.ssym == true ) {
      dev += calc_deviation( *WdT_inactive, *WdT_active );
    }
    W_devstat.add( dev );

#if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : recalculated W "
         << "with deviation = " << dev << endl;

    if ( dev > W_devstat.target ) {
      cout << "HubbardModelVMC::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "HubbardModelVMC::perform_W_update() : approximate W^T =" << endl
           << *WbuT_inactive << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *WdT_inactive << endl;
      }
      cout << "HubbardModelVMC::perform_W_update() : exact W^T =" << endl
           << *WbuT_active << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *WdT_active << endl;
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
    swap( WbuT_inactive, WbuT_active );
    if ( M.ssym == true ) {
      swap( WdT_inactive, WdT_active );
    }

    // updated W should now be in the active buffer
    // debug check recalc W should be in the inactive buffer

    fptype dev = calc_deviation( *WbuT_inactive, *WbuT_active );
    if ( M.ssym == true ) {
      dev += calc_deviation( *WdT_inactive, *WdT_active );
    }

# if VERBOSE >= 1
    cout << "HubbardModelVMC::perform_W_update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > W_devstat.target ) {
      cout << "HubbardModelVMC::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "HubbardModelVMC::perform_W_update() : quickly updated W^T =" << endl
           << *WbuT_active << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *WdT_active << endl;
      }
      cout << "HubbardModelVMC::perform_W_update() : exact W^T =" << endl
           << *WbuT_inactive << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << *WdT_inactive << endl;
      }
    }
# endif
    assert( dev < W_devstat.target );
#endif
  }
}



Eigen::MatrixXfp HubbardModelVMC::calc_DbT() const
{
  assert( M.ssym == false );

  Eigen::MatrixXfp DbT( econf.N(), econf.N() );
  for ( unsigned int eid = 0; eid < econf.N(); ++eid ) {
    DbT.col( eid ) = M.orbitalsT.col( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_DbT() : Db^T = " << endl << DbT << endl;
#endif

  return DbT;
}



Eigen::MatrixXfp HubbardModelVMC::calc_DuT() const
{
  assert( M.ssym == true );

  Eigen::MatrixXfp DuT( econf.N() / 2, econf.N() / 2 );
  for ( unsigned int eid = 0; eid < econf.N() / 2; ++eid ) {
    DuT.col( eid ) = M.orbitalsT.col( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_DuT() : DuT = " << endl << DuT << endl;
#endif

  return DuT;
}



Eigen::MatrixXfp HubbardModelVMC::calc_DdT() const
{
  assert( M.ssym == true );

  Eigen::MatrixXfp DdT( econf.N() / 2, econf.N() / 2 );
  for ( unsigned int eid = econf.N() / 2; eid < econf.N(); ++eid ) {
    DdT.col( eid - econf.N() / 2 )
      = M.orbitalsT.col( lat->get_spinup_site( econf.get_electron_pos( eid ) ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC::calc_DdT() : Dd^T = " << endl << DdT << endl;
#endif

  return DdT;
}


void HubbardModelVMC::calc_new_W()
{
  if ( M.ssym == true ) {

#ifdef USE_ATLAS

    Eigen::MatrixXfp DT = calc_DuT();
    Eigen::Matrix< fptype, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor >
      inMT_outW = M.orbitalsT;
    Eigen::VectorXi ipiv( DT.rows() );
    int info;


    info =
#ifdef USE_FP_DBLPREC
    clapack_dgesv(
#else
    clapack_sgesv(
#endif
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
      CblasRowMajor,
#else
      CblasColMajor,
#endif
      DT.rows(),
      inMT_outW.cols(),
      DT.data(),
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
      DT.cols(),
#else
      DT.rows(),
#endif
      ipiv.data(),
      inMT_outW.data(),
      inMT_outW.rows()
    );
    assert( info == 0 );

    *WbuT_inactive = inMT_outW;


    DT = calc_DdT();
    inMT_outW = M.orbitalsT;

    info =
#ifdef USE_FP_DBLPREC
    clapack_dgesv(
#else
    clapack_sgesv(
#endif
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
      CblasRowMajor,
#else
      CblasColMajor,
#endif
      DT.rows(),
      inMT_outW.cols(),
      DT.data(),
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
      DT.cols(),
#else
      DT.rows(),
#endif
      ipiv.data(),
      inMT_outW.data(),
      inMT_outW.rows()
    );
    assert( info == 0 );

    *WdT_inactive = inMT_outW;

#else

    *WbuT_inactive = calc_DuT().partialPivLu().solve( M.orbitalsT );
    *WdT_inactive  = calc_DdT().partialPivLu().solve( M.orbitalsT );

#endif

    swap( WbuT_inactive, WbuT_active );
    swap( WdT_inactive, WdT_active );

  } else {

#ifdef USE_ATLAS

#else

    *WbuT_inactive = calc_DbT().partialPivLu().solve( M.orbitalsT );

#endif

    swap( WbuT_inactive, WbuT_active );

  }
}



void HubbardModelVMC::calc_qupdated_W( const ElectronHop& hop )
{
  unsigned int k     = hop.k;
  unsigned int l     = hop.l;
  unsigned int k_pos = hop.k_pos;

  Eigen::MatrixXfp*& WT = ( M.ssym == true && hop.k >= econf.N() / 2 ) ?
                          WdT_active : WbuT_active;

  if ( M.ssym == true && hop.k >= econf.N() / 2 ) {
    k     -= econf.N() / 2;
    l     -= lat->L;
    k_pos -= lat->L;
  }

#ifdef USE_ATLAS

  tempWTrow = WT->row( k );
  tempWTcol = WT->col( l ) - WT->col( k_pos );

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
    WT->rows(),
    WT->cols(),
    - 1.f / ( *WT )( k, l ),
    tempWTcol.data(),
    1,
    tempWTrow.data(),
    1,
    WT->data(),
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
    WT->cols()
#else
    WT->rows()
#endif
  );

#else // #ifndef USE_ATLAS

  Eigen::MatrixXfp*& WT_inactive = ( M.ssym == true && hop.k >= econf.N() / 2 ) ?
                                   WdT_inactive : WbuT_inactive;

  *WT_inactive = *WT;

  WT_inactive->noalias() -=
    ( WT->col( l ) - WT->col( k_pos ) )
    * ( WT->row( k ) / ( *WT )( k, l ) );

  swap( WT_inactive, WT );

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
          if ( M.ssym == true && k >= econf.N() / 2 ) {
            sum_Xnn += R_j * ( *WdT_active )( k - econf.N() / 2, *l_it - lat->L );
          } else {
            sum_Xnn += R_j * ( *WbuT_active )( k, *l_it );
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
