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

#include "lattice_2dsquare.hpp"

#include <algorithm>

#if VERBOSE >= 1
# include <iostream>
#endif

# include "utils.hpp"

using namespace std;



Lattice2DSquare::Lattice2DSquare( unsigned int L_init )
  : Lattice( L_init ), S( uintsqrt( L_init ) )
{
  assert( is_perfect_square( L_init ) );
}



int Lattice2DSquare::x( Lattice::spindex l ) const
{
  assert( l < 2 * L );
  // we don't have to convert the spindex to an index first, because
  // it does not make any difference for the position in x
  return l % S;
}

int Lattice2DSquare::y( Lattice::spindex l ) const
{
  assert( l < 2 * L );
  return get_index_from_spindex( l ) / S;
}

int Lattice2DSquare::d( int p1, int p2 ) const
{
  assert( p1 >= 0 && p2 >= 0 );

  int dist = p2 - p1;

  // wrap large d around the boundaries
  if ( dist > static_cast<int>( S ) / 2 ) {
    dist -= S ;
  }
  if ( dist < -1 * static_cast<int>( S ) / 2 ) {
    dist += S;
  }

  return dist;
}



void Lattice2DSquare::get_Xnn(
  Lattice::spindex l, unsigned int X, std::vector<Lattice::spindex>* outbuf ) const
{
  assert( l < 2 * L );
  assert( X == 1 || X == 2 || X == 3 );

  if ( X == 1 ) {
    get_1nn( l, outbuf );
  } else if ( X == 2 ) {
    get_2nn( l, outbuf );
  } else { /* X == 3 */
    assert( X == 3 );
    get_3nn( l, outbuf );
  }
}

void Lattice2DSquare::get_1nn(
  Lattice::spindex l, vector<Lattice::spindex>* outbuf ) const
{
  outbuf->resize( 4 );

  const unsigned int xl = x( l );
  const unsigned int yl = y( l );

  // add left neighbor
  if ( xl == 0 ) {
    ( *outbuf )[0] = l + S - 1;
  } else {
    ( *outbuf )[0] = l - 1;
  }

  // add right neighbor
  if ( xl == S - 1 ) {
    ( *outbuf )[1] = l + 1 - S;
  } else {
    ( *outbuf )[1] = l + 1;
  }

  // add bottom neighbor
  if ( yl == 0 ) {
    ( *outbuf )[2] = l + L - S;
  } else {
    ( *outbuf )[2] = l - S;
  }

  // add top neighbor
  if ( yl == S - 1 ) {
    ( *outbuf )[3] = l + S - L;
  } else {
    ( *outbuf )[3] = l + S;
  }
}

void Lattice2DSquare::get_2nn(
  Lattice::spindex l, vector<Lattice::spindex>* outbuf ) const
{
  outbuf->resize( 4 );

  const unsigned int xl = x( l );
  const unsigned int yl = y( l );

  // add bottom left neighbor
  if ( xl == 0 ) {
    // in left column
    if ( yl == 0 ) {
      // bottom left corner
      ( *outbuf )[0] = l + L - 1;
    } else {
      ( *outbuf )[0] = l - 1;
    }
  } else if ( yl == 0 ) {
    // in bottom row
    // (but NOT in bottom left corner!)
    ( *outbuf )[0] = l + L - S - 1;
  } else {
    // in the center
    ( *outbuf )[0] = l - S - 1;
  }

  // add bottom right neighbor
  if ( xl == S - 1 ) {
    // in right column
    if ( yl == 0 ) {
      // bottom right corner
      ( *outbuf )[1] = l + L + 1 - 2 * S;
    } else {
      ( *outbuf )[1] = l + 1 - 2 * S;
    }
  } else if ( yl == 0 ) {
    // in bottom row
    // (but NOT in bottom right corner!)
    ( *outbuf )[1] = l + L + 1 - S;
  } else {
    // in the center
    ( *outbuf )[1] = l + 1 - S;
  }

  // add top right neighbor
  if ( xl == S - 1 ) {
    // in right column
    if ( yl == S - 1 ) {
      // top right corner
      ( *outbuf )[2] = l + 1 - L;
    } else {
      ( *outbuf )[2] = l + 1;
    }
  } else if ( yl == S - 1 ) {
    // in top row
    // (but NOT in top right corner!)
    ( *outbuf )[2] = l + S + 1 - L;
  } else {
    // in the center
    ( *outbuf )[2] = l + S + 1;
  }

  // add top left neighbor
  if ( xl == 0 ) {
    // in left column
    if ( yl == S - 1 ) {
      // top left corner
      ( *outbuf )[3] = l + 2 * S - 1 - L ;
    } else {
      ( *outbuf )[3] = l + 2 * S - 1;
    }
  } else if ( yl == S - 1 ) {
    // in top row
    // (but NOT in top left corner!)
    ( *outbuf )[3] = l + S - 1 - L;
  } else {
    // in the center
    ( *outbuf )[3] = l + S - 1;
  }
}

void Lattice2DSquare::get_3nn(
  Lattice::spindex l, vector<Lattice::spindex>* outbuf ) const
{
  outbuf->resize( 4 );

  const unsigned int xl = x( l );
  const unsigned int yl = y( l );

  // add left neighbor
  if ( xl <= 1  ) {
    ( *outbuf )[0] = l + S - 2;
  } else {
    ( *outbuf )[0] = l - 2;
  }

  // add right neighbor
  if ( xl >= S - 2 ) {
    ( *outbuf )[1] = l + 2 - S;
  } else {
    ( *outbuf )[1] = l + 2;
  }

  // add bottom neighbor
  if ( yl <= 1 ) {
    ( *outbuf )[2] = l + L - 2 * S;
  } else {
    ( *outbuf )[2] = l - 2 * S;
  }

  // add top neighbor
  if ( yl >= S - 2 ) {
    ( *outbuf )[3] = l + 2 * S - L;
  } else {
    ( *outbuf )[3] = l + 2 * S;
  }
}



Lattice::irridxrel Lattice2DSquare::reduce_idxrel(
  Lattice::spindex i, Lattice::spindex j ) const
{
  assert( i < 2 * L );
  assert( j < 2 * L );
  assert( get_spindex_type( i ) == get_spindex_type( j ) );

  // calculate the absolute values of the position differences
  unsigned int dx = std::abs( d( x( i ), x( j ) ) );
  unsigned int dy = std::abs( d( y( i ), y( j ) ) );

  // dx should be larger than dy
  if ( dy > dx ) {
    swap( dx, dy );
  }

  return dx + S * dy;
}



set<Lattice::irridxrel> Lattice2DSquare::get_all_irridxrels() const
{
  set<Lattice::irridxrel> allrels;
  for ( Lattice::index i = 0; i < L; ++i ) {
    allrels.insert( reduce_idxrel( 0, i ) );
  }

#if VERBOSE >= 1
  cout << "Lattice2DSquare::irreducible_idxrel_list() : "
       << "list of irreducible index relations =" << endl;
  for ( auto it = allrels.begin(); it != allrels.end(); ++it ) {
    cout << *it << endl;
  }
#endif

  return allrels;
}



Lattice::irridxrel Lattice2DSquare::get_maxdist_irridxrel() const
{
  if ( S % 2 == 0 ) {
    return L / 2 + S / 2;
  } else {
    return L / 2;
  }
}



Eigen::VectorXd Lattice2DSquare::r( Lattice::index i, Lattice::index j ) const
{
  assert( i < L );
  assert( j < L );

  Eigen::VectorXd result = Eigen::VectorXd::Zero( 2 );
  result( 0 ) = x( j ) - x( i );
  result( 1 ) = y( j ) - y( i );
  return result;
}



vector<Eigen::VectorXd> Lattice2DSquare::get_qvectors() const
{
  vector<Eigen::VectorXd> allq;
  allq.reserve( S * S / 4 );

  for ( unsigned int i = 0; i <= S / 2; ++i ) {
    for ( unsigned int l = 0; l <= S / 2; ++l ) {
      if ( i == 0 && l == 0 ) {
        continue;
      }
      Eigen::VectorXd q( 2 );
      q[0] = i * 2.0 * M_PI / static_cast<double>( S );
      q[1] = l * 2.0 * M_PI / static_cast<double>( S );
      allq.push_back( q );
    }
  }

  return allq;
}



double Lattice2DSquare::pairsym_modifier(
  optpairsym_t sym, Lattice::spindex i, Lattice::spindex j ) const
{
  assert( get_spindex_type( i ) == get_spindex_type( j ) );

  if ( sym == OPTION_PAIRING_SYMMETRY_SWAVE ) {
    return 1.0;
  } else { // dwave or twisted dwave (doesn't make a difference for 1 plane)

    int dx = d( x( i ), x( j ) );
    int dy = d( y( i ), y( j ) );

    if ( ( dx == 1 || dx == -1 ) && dy == 0 ) {
      // nearest neighbors along x-axis
      return 1.0;
    } else if ( dx == 0 && ( dy == 1 || dy == -1 ) ) {
      // nearest neighbors along y-axis
      return -1.0;
    } else if ( ( dx == 1 || dx == -1 ) && ( dy == 1 || dy == -1 ) ) {
      // second nearest neighbors along diagonal
      return 0.0;
    } else if ( ( dx == 2 || dx == -2 ) && dy == 0 ) {
      // third nearest neighbors along x-axis
      return 1.0;
    } else if ( dx == 0 && ( dy == 2 || dy == -2 ) ) {
      // third nearest neighbors along y-axis
      return -1.0;
    }
  }

  // still here? no meaningful decision yet??? --> this is a bug ...
  assert( false );
  return 0.0; // <-- should never be reached; only to suppress compiler warning
}


unsigned int Lattice2DSquare::get_index_sublattice( Lattice::index i ) const
{
  assert( i < L );
  return ( x( i ) + y( i ) ) % 2;
}
