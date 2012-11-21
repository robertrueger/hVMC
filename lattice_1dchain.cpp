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

#include "lattice_1dchain.hpp"
using namespace std;



Lattice1DChain::Lattice1DChain( unsigned int L_init )
  : Lattice( LATTICE_1DCHAIN, L_init ) { }



vector<unsigned int> Lattice1DChain::get_Xnn( unsigned int l, unsigned int X ) const
{
  assert( l < 2 * L );
  assert( X == 1 || X == 2 || X == 3 );

  vector<unsigned int> Xnn( 2 );

  // add left neighbor to the list
  if ( l < X || ( l >= L && l < L + X ) ) {
    Xnn[0] = l + L - X;
  } else {
    Xnn[0] = l - X;
  }

  // add right neighbor to the list
  if ( ( l < L && l >= L - X ) || l >= 2 * L - X ) {
    Xnn[1] = l + X - L;
  } else {
    Xnn[1] = l + X;
  }

  return Xnn;
}



IrreducibleIdxRel Lattice1DChain::reduce_idxrel(
  unsigned int i, unsigned int j ) const
{
  assert( i < 2 * L );
  assert( j < 2 * L );
  assert( ( i < L && j < L ) || ( i >= L && j >= L ) );

  unsigned int d = j > i ? j - i : i - j;
  IrreducibleIdxRel result( 0, min( d, L - d ) );

#if VERBOSE >= 2
  cout << "Lattice1DChain::reduce_idxrel() : reduction "
       << "(" << i << "," << j << ") -> "
       << "(" << result.first << "," << result.second << ")" << endl;
#endif

  assert( irreducible_idxrel_list().count( result ) == 1 );

  return result;
}



set<IrreducibleIdxRel> Lattice1DChain::irreducible_idxrel_list() const
{
  set<IrreducibleIdxRel> irr_idxrels;
  for ( unsigned int d = 0; d <= L - d; ++d ) {
    assert( d < L );
    irr_idxrels.insert( IrreducibleIdxRel( 0, d ) );
  }

#if VERBOSE >= 2
  cout << "Lattice1DChain::irreducible_idxrel_list() : "
       << "list of irreducible index relations =" << endl;
  for ( auto it = irr_idxrels.begin(); it != irr_idxrels.end(); ++it ) {
    cout << "(" << it->first << "," << it->second << ")" << endl;
  }
#endif

  return irr_idxrels;
}
