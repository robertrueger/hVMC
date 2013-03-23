/*
 * Copyright (c) 2013, Robert Rueger <rueger@itp.uni-frankfurt.de>
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

#include "wmatrix.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <algorithm>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/LU>

#ifdef USE_CBLAS
extern "C" {
# include <cblas.h>
}
#endif

using namespace std;


WMatrix::WMatrix(
  const Lattice* lat_init,
  const SingleParticleOrbitals& detwf_init,
  const ElectronConfiguration& econf_init,
  fptype deviation_target,
  unsigned int updates_until_recalc_init )
  : lat( lat_init ), detwf( detwf_init ), econf( econf_init ),
    Wbu_1(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, econf.N() / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, econf.N() )
    ),
    Wbu_2(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, econf.N() / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, econf.N() )
    ),
    Wbu_active(   &Wbu_1 ),
    Wbu_inactive( &Wbu_2 ),
    Wd_1(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, econf.N() / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    Wd_2(
      detwf_init.ssym ?
      Eigen::MatrixXfp( lat->L, econf.N() / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    Wd_active(   detwf.ssym ? &Wd_1 : nullptr ),
    Wd_inactive( detwf.ssym ? &Wd_2 : nullptr ),
#ifdef USE_CBLAS
    tempWcol(
      detwf.ssym ?
      Eigen::VectorXfp( lat->L ) :
      Eigen::VectorXfp( 2 * lat->L )
    ),
    tempWrow(
      detwf.ssym ?
      Eigen::VectorXfp( econf.N() / 2 ) :
      Eigen::VectorXfp( econf.N() )
    ),
#endif
    updates_until_recalc( updates_until_recalc_init ),
    updates_since_recalc( 0 ),
    devstat( FPDevStat( deviation_target ) ) { }



bool WMatrix::init_and_check()
{
  if ( detwf.ssym == true ) {

    // check spin up part
    Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Du().transpose() );

    if ( lu_decomp.isInvertible() != true ) {
#if VERBOSE >= 2
      cout << "WMatrix::WMatrix() : spin up part has no "
           << "overlap with the determinantal wavefunction" << endl;
#endif
      return false;
    }

    Wbu_active->noalias()
      = lu_decomp.solve( detwf.M.transpose() ).transpose();
    fptype Wbu_avg
      = Wbu_active->squaredNorm() / static_cast<fptype>( Wbu_active->size() );
    if ( Wbu_avg > 100.f ) {
#if VERBOSE >= 2
      cout << "WMatrix::WMatrix() : spin up part has too "
           << "little overlap with the determinantal wavefunction, "
           << "inverse overlap measure is: " << Wbu_avg << endl;
#endif
      return false;
    }

    // check spin down part
    lu_decomp.compute( calc_Dd().transpose() );
    if ( lu_decomp.isInvertible() != true ) {
#if VERBOSE >= 2
      cout << "WMatrix::WMatrix() : spin down part has no "
           << "overlap with the determinantal wavefunction" << endl;
#endif
      return false;
    }

    Wd_active->noalias()
      = lu_decomp.solve( detwf.M.transpose() ).transpose();
    fptype Wd_avg
      = Wd_active->squaredNorm()  / static_cast<fptype>( Wd_active->size() );
    if ( Wd_avg > 100.f ) {
#if VERBOSE >= 2
      cout << "WMatrix::WMatrix() : spin down part has too "
           << "little overlap with the determinantal wavefunction, "
           << "inverse overlap measure is: " << Wd_avg << endl;
#endif
      return false;
    }

  } else {

    // check whole determinantal part
    Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Db().transpose() );
    if ( lu_decomp.isInvertible() != true ) {
#if VERBOSE >= 2
      cout << "WMatrix::WMatrix() : state has no "
           << "overlap with the determinantal wavefunction" << endl;
#endif
      return false;
    }

    Wbu_active->noalias()
      = lu_decomp.solve( detwf.M.transpose() ).transpose();
    fptype Wbu_avg
      = Wbu_active->squaredNorm() / static_cast<fptype>( Wbu_active->size() );
    if ( Wbu_avg > 50.f ) {
#if VERBOSE >= 2
      cout << "WMatrix::WMatrix() : state has too "
           << "little overlap with the determinantal wavefunction, "
           << "inverse overlap measure is: " << Wbu_avg << endl;
#endif
      return false;
    }

  }

  return true;
}



fptype WMatrix::operator()( unsigned int i, unsigned int j ) const
{
  assert( i < 2 * lat->L && j < econf.N() );

  if ( detwf.ssym == true ) {
    assert(
      ( i < lat->L && j < econf.N() / 2 ) ||
      ( i >= lat->L && j >= econf.N() / 2 )
    );

    if ( i >= lat->L ) {
      return ( *Wd_active  )( i - lat->L, j - econf.N() / 2 );
    } else {
      return ( *Wbu_active )( i, j );
    }

  } else {

    return ( *Wbu_active )( i, j );

  }
}



void WMatrix::update( const ElectronHop& hop )
{
  if ( updates_since_recalc >= updates_until_recalc ) {

#if VERBOSE >= 2
    cout << "WMatrix::perform_W_update() : recalculating W!" << endl;
#endif

    updates_since_recalc = 0;

    // puts updated W into the active buffer
    // (only buffer of the hopping spin direction is changed)
    calc_qupdated( hop );

    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new();

    fptype dev = calc_deviation( *Wbu_inactive, *Wbu_active );
    if ( detwf.ssym == true ) {
      dev += calc_deviation( *Wd_inactive, *Wd_active );
    }
    devstat.add( dev );

#if VERBOSE >= 2
    cout << "WMatrix::perform_W_update() : recalculated W "
         << "with deviation = " << dev << endl;

    if ( dev > devstat.target ) {
      cout << "WMatrix::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "WMatrix::perform_W_update() : approximate W =" << endl
           << *Wbu_inactive << endl;
      if ( detwf.ssym == true ) {
        cout << "----->" << endl << *Wd_inactive << endl;
      }
      cout << "WMatrix::perform_W_update() : exact W =" << endl
           << *Wbu_active << endl;
      if ( detwf.ssym == true ) {
        cout << "----->" << endl << *Wd_active << endl;
      }
    }
#endif

    assert( dev < devstat.target );

  } else {

#if VERBOSE >= 2
    cout << "WMatrix::perform_W_update() : "
         << "performing a quick update of W!" << endl;
#endif

    ++updates_since_recalc;

    // puts updated W into the active buffer
    // (only buffer of the hopping spin direction is changed)
    calc_qupdated( hop );

#ifndef NDEBUG
    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new();

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

# if VERBOSE >= 2
    cout << "WMatrix::perform_W_update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > devstat.target ) {
      cout << "WMatrix::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "WMatrix::perform_W_update() : quickly updated W =" << endl
           << *Wbu_active << endl;
      if ( detwf.ssym == true ) {
        cout << "----->" << endl << *Wd_active << endl;
      }
      cout << "WMatrix::perform_W_update() : exact W =" << endl
           << *Wbu_inactive << endl;
      if ( detwf.ssym == true ) {
        cout << "----->" << endl << *Wd_inactive << endl;
      }
    }
# endif
    assert( dev < devstat.target );
#endif
  }
}



void WMatrix::calc_new()
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



void WMatrix::calc_qupdated( const ElectronHop& hop )
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

  Eigen::MatrixXfp*& W_inactive
    = ( detwf.ssym == true && hop.k >= econf.N() / 2 ) ?
      Wd_inactive : Wbu_inactive;

  *W_inactive = *W;

  W_inactive->noalias() -=
    ( W->col( k ) / ( *W )( l, k ) )
    * ( W->row( l ) - W->row( k_pos ) );

  swap( W_inactive, W );

#endif
}



Eigen::MatrixXfp WMatrix::calc_Db() const
{
  assert( detwf.ssym == false );

  Eigen::MatrixXfp Db( econf.N(), econf.N() );
  for ( unsigned int eid = 0; eid < econf.N(); ++eid ) {
    Db.row( eid ) = detwf.M.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 3
  cout << "WMatrix::calc_Db() : Db = " << endl << Db << endl;
#endif

  return Db;
}



Eigen::MatrixXfp WMatrix::calc_Du() const
{
  assert( detwf.ssym == true );

  Eigen::MatrixXfp Du( econf.N() / 2, econf.N() / 2 );
  for ( unsigned int eid = 0; eid < econf.N() / 2; ++eid ) {
    Du.row( eid ) = detwf.M.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 3
  cout << "WMatrix::calc_Du() : Du = " << endl << Du << endl;
#endif

  return Du;
}



Eigen::MatrixXfp WMatrix::calc_Dd() const
{
  assert( detwf.ssym == true );

  Eigen::MatrixXfp Dd( econf.N() / 2, econf.N() / 2 );
  for ( unsigned int eid = econf.N() / 2; eid < econf.N(); ++eid ) {
    Dd.row( eid - econf.N() / 2 )
      = detwf.M.row( lat->get_spinup_site( econf.get_electron_pos( eid ) ) );
  }

#if VERBOSE >= 3
  cout << "WMatrix::calc_Dd() : Dd = " << endl << Dd << endl;
#endif

  return Dd;
}



FPDevStat WMatrix::get_devstat() const
{
  return devstat;
}
