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

#if VERBOSE >= 2
# include <iostream>
#endif

#include <cmath>
#include <algorithm>

#include "macros.h"

using namespace std;



TVector::TVector(
  const Lattice* lat_init,
  const Jastrow& v_init,
  const ParticleConfiguration& pconf_init,
  double deviation_target,
  unsigned int updates_until_recalc_init )
  : lat( lat_init ), v( v_init ), pconf( pconf_init ),
    T( lat->L ),
    updates_until_recalc( updates_until_recalc_init ),
    updates_since_recalc( 0 ),
    devstat( FPDevStat( deviation_target ) ) { }



void TVector::init()
{
  T = calc_new();
  devstat.reset();
  updates_since_recalc = 0;

#if VERBOSE >= 2
    cout << "TVector::init() : initial T = " << endl << T.transpose() << endl;
#endif
}



const Eigen::VectorXd& TVector::get() const
{
  return T;
}



void TVector::update( const ParticleHop& hop )
{
  if ( updates_since_recalc >= updates_until_recalc ) {

#if VERBOSE >= 2
    cout << "TVector::update() : recalculating T!" << endl;
#endif

    updates_since_recalc = 0;

    Eigen::VectorXd T_approx = calc_qupdated( hop );
    T = calc_new();

    double dev = calc_deviation( T_approx, T );
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
    Eigen::VectorXd T_chk = calc_new();

    double dev = calc_deviation( T, T_chk );

# if VERBOSE >= 2
    cout << "TVector::update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > devstat.target || !std::isfinite( dev ) ) {
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



Eigen::VectorXd TVector::calc_new() const
{
  Eigen::VectorXd T_new = Eigen::VectorXd::Zero( lat->L );

  for ( Lattice::index i = 0; i < lat->L; ++i ) {
    for ( Lattice::index j = 0; j < lat->L; ++j ) {
      T_new( i ) +=
        v( i, j ) * ( pconf.get_spindex_occ()[ j ] -
                      pconf.get_spindex_occ()[ j + lat->L ] );
    }
  }
  return T_new;
}



Eigen::VectorXd TVector::calc_qupdated( const ParticleHop& hop ) const
{
  Eigen::VectorXd T_diff( lat->L );

  for ( Lattice::index i = 0; i < lat->L; ++i ) {
    T_diff( i ) = v( i, lat->get_index_from_spindex( hop.l ) )
                  - v( i, lat->get_index_from_spindex( hop.k ) );
  }

  return T + ( hop.l < lat->L ? 1.0 : -1.0 ) * T_diff;
}



FPDevStat TVector::get_devstat() const
{
  return devstat;
}
