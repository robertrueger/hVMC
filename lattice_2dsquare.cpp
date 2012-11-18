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

#include "lattice_2dsquare.hpp"
using namespace std;



Lattice2DSquare::Lattice2DSquare( unsigned int L_init )
  : Lattice( LATTICE_2DSQUARE, L_init ), S( uintsqrt( L_init ) )
{
  assert( is_perfect_square( L_init ) );
}



vector<unsigned int> Lattice2DSquare::get_Xnn( unsigned int l, unsigned int X ) const
{
  assert( l < 2 * L );
  assert( X == 1 || X == 2 || X == 3 );

  // TODO: optimize by caching position bools (Robert Rueger, 2012-11-18 00:53)

  if ( X == 1 ) {
    return get_1nn( l, X );
  } else if ( X == 2 ) {
    return get_2nn( l, X );
  } else /* X == 3 */ {
    return get_3nn( l, X );
  }
}

vector<unsigned int> Lattice2DSquare::get_1nn( unsigned int l, unsigned int X ) const
{
  vector<unsigned int> nn1( 4 );

  // add left neighbor
  if ( l % S == 0 ) {
    nn1[0] = l + S - 1;
  } else {
    nn1[0] = l - 1;
  }

  // add right neighbor
  if ( ( l + 1 ) % S == 0 ) {
    nn1[1] = l + 1 - S;
  } else {
    nn1[1] = l + 1;
  }

  // add bottom neighbor
  if ( l < S || ( l >= L && l < L + S ) ) {
    nn1[2] = l + L - S;
  } else {
    nn1[2] = l - S;
  }

  // add top neighbor
  if ( l >= 2 * L - S || ( l < L && l >= L - S ) ) {
    nn1[3] = l + S - L;
  } else {
    nn1[3] = l + S;
  }

  return nn1;
}

vector<unsigned int> Lattice2DSquare::get_2nn( unsigned int l, unsigned int X ) const
{
  vector<unsigned int> nn2( 4 );

  // add bottom left neighbor
  if ( l % S == 0 ) {
    // in left column
    if ( l == 0 || l == L ) {
      // bottom left corner
      nn2[0] = l + L - 1;
    } else {
      nn2[0] = l - 1;
    }
  } else if ( l < S || ( l >= L && l < L + S ) ) {
    // in bottom row
    // (but NOT in bottom left corner!)
    nn2[0] = l + L - S - 1;
  } else {
    // in the center
    nn2[0] = l - S - 1;
  }

  // add bottom right neighbor
  if ( ( l + 1 ) % S == 0 ) {
    // in right column
    if ( l == S - 1 || l == L + S - 1 ) {
      // bottom right corner
      nn2[1] = l + L + 1 - 2 * S;
    } else {
      nn2[1] = l + 1 - 2 * S;
    }
  } else if ( l < S || ( l >= L && l < L + S ) ) {
    // in bottom row
    // (but NOT in bottom right corner!)
    nn2[1] = l + L + 1 - S;
  } else {
    // in the center
    nn2[1] = l + 1 - S;
  }

  // add top right neighbor
  if ( ( l + 1 ) % S == 0 ) {
    // in right column
    if ( l == L - 1 || l == 2 * L - 1 ) {
      // top right corner
      nn2[2] = l + 1 - L;
    } else {
      nn2[2] = l + 1;
    }
  } else if ( l >= 2 * L - S || ( l < L && l >= L - S ) ) {
    // in top row
    // (but NOT in top right corner!)
    nn2[2] = l + S + 1 - L;
  } else {
    // in the center
    nn2[2] = l + S + 1;
  }

  // add top left neighbor
  if ( l % S == 0 ) {
    // in left column
    if ( l == L - S || l == 2 * L - S ) {
      // top left corner
      nn2[3] = l + 2 * S - 1 - L ;
    } else {
      nn2[3] = l + 2 * S - 1;
    }
  } else if ( l >= 2 * L - S || ( l < L && l >= L - S ) ) {
    // in top row
    // (but NOT in top left corner!)
    nn2[3] = l + S - 1 - L;
  } else {
    // in the center
    nn2[3] = l + S - 1;
  }

  return nn2;
}

vector<unsigned int> Lattice2DSquare::get_3nn( unsigned int l, unsigned int X ) const
{
  vector<unsigned int> nn3( 4 );

  // add left neighbor
  if ( l % S <= 1  ) {
    nn3[0] = l + S - 2;
  } else {
    nn3[0] = l - 2;
  }

  // add right neighbor
  if ( l % S >= S - 2 ) {
    nn3[1] = l + 2 - S;
  } else {
    nn3[1] = l + 2;
  }

  // add bottom neighbor
  if ( l < 2 * S || ( l >= L && l < L + 2 * S ) ) {
    nn3[2] = l + L - 2 * S;
  } else {
    nn3[2] = l - 2 * S;
  }

  // add top neighbor
  if ( l >= 2 * ( L - S ) || ( l < L && l >= L - 2 * S ) ) {
    nn3[3] = l + 2 * S - L;
  } else {
    nn3[3] = l + 2 * S;
  }

  return nn3;
}



IrreducibleIdxRel Lattice2DSquare::reduce_idxrel(
  unsigned int i, unsigned int j ) const
{
  assert( i < 2 * L );
  assert( j < 2 * L );
  assert( ( i < L && j < L ) || ( i >= L && j >= L ) );


}



set<IrreducibleIdxRel> Lattice2DSquare::irreducible_idxrel_list() const
{
}
