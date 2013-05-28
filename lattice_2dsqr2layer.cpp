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

#include "lattice_2dsqr2layer.hpp"

#include <algorithm>

#if VERBOSE >= 1
# include <iostream>
#endif

#include "utils.hpp"

using namespace std;



Lattice2DSquare2Layer::Lattice2DSquare2Layer( unsigned int L_init )
  : Lattice( L_init ), L_layer( L_init / 2 ), S( uintsqrt( L_init / 2 ) )
{
  assert( L_init % 2 == 0 );
  assert( is_perfect_square( L_init / 2 ) );
}



int Lattice2DSquare2Layer::x( Lattice::spindex l ) const
{
  assert( l < 2 * L );
  // we don't have to convert the spindex to an index first
  // or care about which plane we are in, because
  // it does not make any difference for the position in x
  return l % S;
}

int Lattice2DSquare2Layer::y( Lattice::spindex l ) const
{
  assert( l < 2 * L );
  return ( l % L_layer ) / S;
}

int Lattice2DSquare2Layer::z( Lattice::spindex l ) const
{
  assert( l < 2 * L );
  return ( l % L ) / L_layer;
}

int Lattice2DSquare2Layer::d( int p1, int p2 ) const
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



void Lattice2DSquare2Layer::get_Xnn(
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
    // corresponding sites in the other plane are treated as 3rd n.n.
    get_opn( l, outbuf );
  }
}

void Lattice2DSquare2Layer::get_1nn(
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
    ( *outbuf )[2] = l + L_layer - S;
  } else {
    ( *outbuf )[2] = l - S;
  }

  // add top neighbor
  if ( yl == S - 1 ) {
    ( *outbuf )[3] = l + S - L_layer;
  } else {
    ( *outbuf )[3] = l + S;
  }
}

void Lattice2DSquare2Layer::get_2nn(
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
      ( *outbuf )[0] = l + L_layer - 1;
    } else {
      ( *outbuf )[0] = l - 1;
    }
  } else if ( yl == 0 ) {
    // in bottom row
    // (but NOT in bottom left corner!)
    ( *outbuf )[0] = l + L_layer - S - 1;
  } else {
    // in the center
    ( *outbuf )[0] = l - S - 1;
  }

  // add bottom right neighbor
  if ( xl == S - 1 ) {
    // in right column
    if ( yl == 0 ) {
      // bottom right corner
      ( *outbuf )[1] = l + L_layer + 1 - 2 * S;
    } else {
      ( *outbuf )[1] = l + 1 - 2 * S;
    }
  } else if ( yl == 0 ) {
    // in bottom row
    // (but NOT in bottom right corner!)
    ( *outbuf )[1] = l + L_layer + 1 - S;
  } else {
    // in the center
    ( *outbuf )[1] = l + 1 - S;
  }

  // add top right neighbor
  if ( xl == S - 1 ) {
    // in right column
    if ( yl == S - 1 ) {
      // top right corner
      ( *outbuf )[2] = l + 1 - L_layer;
    } else {
      ( *outbuf )[2] = l + 1;
    }
  } else if ( yl == S - 1 ) {
    // in top row
    // (but NOT in top right corner!)
    ( *outbuf )[2] = l + S + 1 - L_layer;
  } else {
    // in the center
    ( *outbuf )[2] = l + S + 1;
  }

  // add top left neighbor
  if ( xl == 0 ) {
    // in left column
    if ( yl == S - 1 ) {
      // top left corner
      ( *outbuf )[3] = l + 2 * S - 1 - L_layer;
    } else {
      ( *outbuf )[3] = l + 2 * S - 1;
    }
  } else if ( yl == S - 1 ) {
    // in top row
    // (but NOT in top left corner!)
    ( *outbuf )[3] = l + S - 1 - L_layer;
  } else {
    // in the center
    ( *outbuf )[3] = l + S - 1;
  }
}

void Lattice2DSquare2Layer::get_opn(
  Lattice::spindex l, vector<Lattice::spindex>* outbuf ) const
{
  outbuf->resize( 1 );

  if ( z( l ) == 0 ) {
    ( *outbuf )[0] = l + L_layer;
  } else {
    assert( z( l ) == 1 );
    ( *outbuf )[0] = l - L_layer;
  }
}



Lattice::irridxrel Lattice2DSquare2Layer::reduce_idxrel(
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

  // calculate the irreducible index relation if i and j were in the same plane
  Lattice::irridxrel iir = dx + S * dy;

  // if they are in different planes, we obviously need a different
  // irreducible index relation as if they were in the same plane
  if ( z( i ) != z( j ) ) {
    iir += L_layer;
  }

  return iir;
}



set<Lattice::irridxrel> Lattice2DSquare2Layer::get_all_irridxrels() const
{
  set<Lattice::irridxrel> allrels;
  for ( Lattice::index i = 0; i < L; ++i ) {
    allrels.insert( reduce_idxrel( 0, i ) );
  }

#if VERBOSE >= 1
  cout << "Lattice2DSquare2Layer::irreducible_idxrel_list() : "
       << "list of irreducible index relations =" << endl;
  for ( auto it = allrels.begin(); it != allrels.end(); ++it ) {
    cout << *it << endl;
  }
#endif

  return allrels;
}



Lattice::irridxrel Lattice2DSquare2Layer::get_maxdist_irridxrel() const
{
  if ( S % 2 == 0 ) {
    return L_layer / 2 + S / 2 /*+ L_layer*/;
  } else {
    return L_layer / 2 /*+ L_layer*/;
  }
}



Eigen::VectorXd Lattice2DSquare2Layer::r( Lattice::index i, Lattice::index j ) const
{
  assert( i < L );
  assert( j < L );

  Eigen::VectorXd result = Eigen::VectorXd::Zero( 2 );
  result( 0 ) = x( j ) - x( i );
  result( 1 ) = y( j ) - y( i );
  return result;
}



vector<Eigen::VectorXd> Lattice2DSquare2Layer::get_qvectors() const
{
  vector<Eigen::VectorXd> allq;
  allq.reserve( S * S / 4 );

  for ( unsigned int i = 1; i <= S / 2; ++i ) {
    for ( unsigned int l = 1; l <= S / 2; ++l ) {
      Eigen::VectorXd q( 2 );
      q[0] = i * 2.0 * M_PI / static_cast<double>( S );
      q[1] = l * 2.0 * M_PI / static_cast<double>( S );
      allq.push_back( q );
    }
  }

  return allq;
}



double Lattice2DSquare2Layer::pairsym_modifier(
  optpairsym_t sym, Lattice::spindex i, Lattice::spindex j ) const
{
  assert( get_spindex_type( i ) == get_spindex_type( j ) );

  if ( sym == OPTION_PAIRING_SYMMETRY_SWAVE ) {
    return 1.0;
  } else { // dwave or twisted dwave

    int dx = d( x( i ), x( j ) );
    int dy = d( y( i ), y( j ) );

    if ( z( i ) != z( j ) ) {
      assert( dx == 0 && dy == 0 );
      // different planes

      return 1.0;

    } else {
      assert( !( dx == 0 && dy == 0 ) );
      // same plane

      double modifier;

      if ( ( dx == 1 || dx == -1 ) && dy == 0 ) {
        // nearest neighbors along x-axis
        modifier = 1.0;
      } else if ( dx == 0 && ( dy == 1 || dy == -1 ) ) {
        // nearest neighbors along y-axis
        modifier = -1.0;
      } else {
        assert( ( dx == 1 || dx == -1 ) && ( dy == 1 || dy == -1 ) );
        // second nearest neighbors along diagonal
        modifier = 0.0;
      }

      // twisted dwave symmetry reverses the sign in the second plane
      if ( sym == OPTION_PAIRING_SYMMETRY_DWAVE_TWISTED && z( i ) == 1 ) {
        modifier *= -1.0;
      }

      return modifier;
    }
  }
}



Eigen::VectorXi Lattice2DSquare2Layer::get_random_spindex_occ(
   unsigned int Npu, unsigned int Npd, mt19937& rng ) const
{
  // make sure Npu and Npd are even numbers,
  // so that we can distribute them evenly among the planes
  assert( Npu % 2 == 0 && Npd % 2 == 0 );

  Eigen::VectorXi spindex_occ = Eigen::VectorXi::Zero( 2 * L );

  // distribute half of the Npu particles randomly in each plane
  while ( spindex_occ.segment( 0, L_layer ).sum()
            < static_cast<int>( Npu / 2 ) ) {
    spindex_occ[
      uniform_int_distribution<Lattice::index>( 0, L_layer - 1 )( rng )
    ] = 1;
  }
  while ( spindex_occ.segment( L_layer, L_layer ).sum()
            < static_cast<int>( Npu / 2 ) ) {
    spindex_occ[
      uniform_int_distribution<Lattice::index>( L_layer, L - 1 )( rng )
    ] = 1;
  }

  // distribute half of the Npd particles randomly in each plane
  while ( spindex_occ.segment( L, L_layer ).sum()
            < static_cast<int>( Npd / 2 ) ) {
    spindex_occ[
      uniform_int_distribution<Lattice::index>( L, L + L_layer - 1 )( rng )
    ] = 1;
  }
  while ( spindex_occ.segment( L + L_layer, L_layer ).sum()
            < static_cast<int>( Npd / 2 ) ) {
    spindex_occ[
      uniform_int_distribution<Lattice::index>( L + L_layer, 2 * L - 1 )( rng )
    ] = 1;
  }

  return spindex_occ;
}
