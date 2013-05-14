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

#include "macros.h"

using namespace std;



WMatrix::WMatrix(
  const Lattice* lat_init,
  const DeterminantalWavefunction& detwf_init,
  const ParticleConfiguration& pconf_init,
  double deviation_target,
  unsigned int updates_until_recalc_init )
  : lat( lat_init ), detwf( detwf_init ), pconf( pconf_init ),
    W_1( 2 * lat->L, pconf.Np ),
    W_2( 2 * lat->L, pconf.Np ),
    W_active(   &W_1 ),
    W_inactive( &W_2 ),
#ifdef USE_CBLAS
    tempWcol( 2 * lat->L ),
    tempWrow( pconf.Np ),
#endif
    updates_until_recalc( updates_until_recalc_init ),
    updates_since_recalc( 0 ),
    devstat( FPDevStat( deviation_target ) ) { }



bool WMatrix::init_and_check()
{
  Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_D().transpose() );
  if ( lu_decomp.isInvertible() != true ) {
#if VERBOSE >= 1
    cout << "WMatrix::check_and_init() : state has no "
         << "overlap with the determinantal wavefunction, "
         << "D =" << endl << calc_D() << endl;
#endif
    return false;
  }

  W_active->noalias()
    = lu_decomp.solve( detwf.M().transpose() ).transpose();

  return true;
}



fptype WMatrix::operator()( unsigned int i, unsigned int j ) const
{
  assert( i < 2 * lat->L && j < pconf.Np );
  return ( *W_active )( i, j );
}



void WMatrix::update( const ParticleHop& hop )
{
  if ( updates_since_recalc >= updates_until_recalc ) {

#if VERBOSE >= 2
    cout << "WMatrix::update() : recalculating W!" << endl;
#endif

    updates_since_recalc = 0;

    // puts updated W into the active buffer
    // (only buffer of the hopping spin direction is changed)
    calc_qupdated( hop );

    // puts recalculated W in the active buffers
    // (pushs updated W into the inactive buffers)
    calc_new();

    double dev = calc_deviation( *W_inactive, *W_active );
    devstat.add( dev );

#if VERBOSE >= 2
    cout << "WMatrix::update() : recalculated W "
         << "with deviation = " << dev << endl;

    if ( dev > devstat.target ) {
      cout << "WMatrix::update() : deviation goal for matrix "
           << "W not met!" << endl
           << "WMatrix::update() : approximate W =" << endl
           << *W_inactive << endl;
      cout << "WMatrix::update() : exact W =" << endl
           << *W_active << endl;
    }
#endif

    assert( dev < devstat.target );

  } else {

#if VERBOSE >= 2
    cout << "WMatrix::update() : "
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
    swap( W_inactive, W_active );
    // updated W should now be in the active buffer

    // debug check recalculated W
    double dev = calc_deviation( *W_inactive, *W_active );

# if VERBOSE >= 2
    cout << "WMatrix::update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > devstat.target ) {
      cout << "WMatrix::update() : deviation goal for matrix "
           << "W not met!" << endl
           << "WMatrix::update() : quickly updated W =" << endl
           << *W_active << endl;
      cout << "WMatrix::update() : exact W =" << endl
           << *W_inactive << endl;
    }
# endif
    assert( dev < devstat.target );
#endif
  }
}



void WMatrix::calc_new()
{
  W_inactive->noalias()
    = calc_D().transpose().partialPivLu().solve( detwf.M().transpose() ).transpose();

  swap( W_inactive, W_active );
}



void WMatrix::calc_qupdated( const ParticleHop& hop )
{
  unsigned int k         = hop.k;
  Lattice::spindex l     = hop.l;
  Lattice::spindex k_pos = hop.k_pos;

#ifdef USE_CBLAS

  tempWcol = W_active->col( k );
  tempWrow = W_active->row( l ) - W_active->row( k_pos );

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
    W_active->rows(),
    W_active->cols(),
    - 1.f / ( *W_active )( l, k ),
    tempWcol.data(),
    1,
    tempWrow.data(),
    1,
    W_active->data(),
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
    W_active->cols()
#else
    W_active->rows()
#endif
  );

#else // #ifndef USE_CBLAS

  *W_inactive = *W_active;

  W_inactive->noalias() -=
    ( W_active->col( k ) / ( *W_active )( l, k ) )
    * ( W_active->row( l ) - W_active->row( k_pos ) );

  swap( W_inactive, W_active );

#endif
}



Eigen::MatrixXfp WMatrix::calc_D() const
{
  Eigen::MatrixXfp D( pconf.Np, pconf.Np );
  for ( unsigned int pid = 0; pid < pconf.Np; ++pid ) {
    D.row( pid ) = detwf.M().row( pconf.get_particle_pos( pid ) );
  }

#if VERBOSE >= 3
  cout << "WMatrix::calc_D() : D = " << endl << Db << endl;
#endif

  return D;
}



FPDevStat WMatrix::get_devstat() const
{
  return devstat;
}
