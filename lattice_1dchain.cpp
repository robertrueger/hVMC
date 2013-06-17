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

#include "lattice_1dchain.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <cmath>

using namespace std;



Lattice1DChain::Lattice1DChain( unsigned int L_init )
  : Lattice( L_init ) { }



void Lattice1DChain::get_Xnn(
  Lattice::spindex l, unsigned int X, std::vector<Lattice::spindex>* outbuf ) const
{
  assert( l < 2 * L );
  assert( X == 1 || X == 2 || X == 3 );

  outbuf->resize( 2 );

  // add left neighbor to the list
  if ( l < X || ( l >= L && l < L + X ) ) {
    ( *outbuf )[0] = l + L - X;
  } else {
    ( *outbuf )[0] = l - X;
  }

  // add right neighbor to the list
  if ( ( l < L && l >= L - X ) || l >= 2 * L - X ) {
    ( *outbuf )[1] = l + X - L;
  } else {
    ( *outbuf )[1] = l + X;
  }
}



Lattice::irridxrel Lattice1DChain::reduce_idxrel(
  Lattice::spindex i, Lattice::spindex j ) const
{
  assert( i < 2 * L );
  assert( j < 2 * L );
  assert( get_spindex_type( i ) == get_spindex_type( j ) );

  const unsigned int d = j > i ? j - i : i - j;
  const Lattice::irridxrel result = min( d, L - d );

#if VERBOSE >= 2
  cout << "Lattice1DChain::reduce_idxrel() : reduction "
       << "(" << i << "," << j << ") -> " << result << endl;
#endif

  assert( get_all_irridxrels().count( result ) == 1 );

  return result;
}



set<Lattice::irridxrel> Lattice1DChain::get_all_irridxrels() const
{
  set<Lattice::irridxrel> allrels;
  for ( unsigned int d = 0; d <= L - d; ++d ) {
    assert( d < L );
    allrels.insert( d );
  }

#if VERBOSE >= 2
  cout << "Lattice1DChain::irreducible_idxrel_list() : "
       << "list of irreducible index relations =" << endl;
  for ( auto it = allrels.begin(); it != allrels.end(); ++it ) {
    cout << *it << endl;
  }
#endif

  return allrels;
}



Lattice::irridxrel Lattice1DChain::get_maxdist_irridxrel() const
{
  return L / 2;
}



Eigen::VectorXd Lattice1DChain::r( Lattice::index i, Lattice::index j ) const
{
  assert( i < L );
  assert( j < L );

  Eigen::VectorXd result = Eigen::VectorXd::Zero( 1 );
  result( 0 ) = static_cast<double>( j ) - static_cast<double>( i );
  return result;
}



vector<Eigen::VectorXd> Lattice1DChain::get_qvectors() const
{
  vector<Eigen::VectorXd> allq;
  allq.reserve( L / 2 );

  for ( unsigned int i = 1; i <= L / 2; ++i ) {
    Eigen::VectorXd q( 1 );
    q[0] = i * 2.0 * M_PI / static_cast<double>( L );
    allq.push_back( q );
  }

  return allq;
}



double Lattice1DChain::pairsym_modifier(
  optpairsym_t, Lattice::spindex, Lattice::spindex ) const
{
  // s- or d-wave symmetry doesn't make a difference in 1D ...
  return 1.0;
}



unsigned int Lattice1DChain::get_index_sublattice( Lattice::index i ) const
{
  assert( i < L );
  return i % 2;
}
