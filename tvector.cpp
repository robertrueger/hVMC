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

#include "tvector.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <algorithm>

#include "macros.h"

using namespace std;



TVector::TVector(
  const Lattice* lat_init,
  const Jastrow& v_init,
  const ElectronConfiguration& econf_init,
  fptype deviation_target,
  unsigned int updates_until_recalc_init )
  : lat( lat_init ), v( v_init ), econf( econf_init ),
    T( Eigen::VectorXfp( lat->L ) ),
    updates_until_recalc( updates_until_recalc_init ),
    updates_since_recalc( 0 ),
    devstat( FPDevStat( deviation_target ) ) { }



bool TVector::init_and_check()
{
  T = calc_new();
  fptype T_avg = T.squaredNorm() / static_cast<fptype>( T.size() );
  if ( T_avg > 100.f ) {
#if VERBOSE >= 2
    cout << "TMatrix::init_and_check() : Jastrow ratios "
         << "are to small, inverse measure is: " << T_avg << endl;
#endif
    return false;
  } else {
    return true;
  }
}



fptype TVector::operator()( unsigned int i ) const
{
  assert( i < lat->L );
  return T( i );
}



void TVector::update( const ElectronHop& hop )
{
  if ( updates_since_recalc >= updates_until_recalc ) {

#if VERBOSE >= 2
    cout << "TVector::update() : recalculating T!" << endl;
#endif

    updates_since_recalc = 0;

    const Eigen::MatrixXfp& T_approx = calc_qupdated( hop );
    T = calc_new();

    fptype dev = calc_deviation( T_approx, T );
    devstat.add( dev );

#if VERBOSE >= 2
    cout << "TVector::update() : recalculated T "
         << "with deviation = " << dev << endl;

    if ( dev > devstat.target ) {
      cout << "TVector::update() : deviation goal for matrix "
           << "T not met!" << endl
           << "TVector::update() : approximate T =" << endl
           << T_approx.transpose() << endl
           << "TVector::update() : exact T =" << endl
           << T.transpose() << endl;
    }
#endif

    assert( dev < devstat.target );

  } else {

#if VERBOSE >= 2
    cout << "TVector::update() : "
         << "performing a quick update of T!" << endl;
#endif

    ++updates_since_recalc;

    T = calc_qupdated( hop );

#ifndef NDEBUG
    const Eigen::MatrixXfp& T_chk = calc_new();
    fptype dev = calc_deviation( T, T_chk );

# if VERBOSE >= 2
    cout << "TVector::update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > devstat.target ) {
      cout << "TVector::update() : deviation goal for matrix "
           << "T not met!" << endl
           << "TVector::update() : quickly updated T =" << endl
           << T.transpose() << endl
           << "TVector::update() : exact T =" << endl
           << T_chk.transpose() << endl;
    }
# endif
#endif

    assert( dev < devstat.target );
  }
}



Eigen::VectorXfp TVector::calc_new() const
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



Eigen::VectorXfp TVector::calc_qupdated( const ElectronHop& hop ) const
{
  Eigen::VectorXfp T_prime( lat->L );

  for ( unsigned int i = 0; i < lat->L; ++i ) {
    T_prime( i ) = T( i ) * v.exp( i, lat->get_spinup_site( hop.l ) )
                   / v.exp( i, lat->get_spinup_site( hop.k_pos ) );
  }

  return T_prime;
}



FPDevStat TVector::get_devstat() const
{
  return devstat;
}
