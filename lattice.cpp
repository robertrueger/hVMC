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

#include "lattice.hpp"
using namespace std;


Lattice::index Lattice::get_index_from_spindex( Lattice::spindex l ) const
{
  assert( l < 2 * L );

  if ( get_spindex_type( l ) == Lattice::spindex_type::down ) {
    return l - L;
  } else {
    return l;
  }
}


Lattice::spindex Lattice::get_linked_spindex( Lattice::spindex l ) const
{
  assert( l < 2 * L );

  if ( get_spindex_type( l ) == Lattice::spindex_type::up ) {
    return l + L;
  } else {
    return l - L;
  }
}


Lattice::spindex_type Lattice::get_spindex_type( Lattice::spindex l ) const
{
  assert( l < 2 * L );

  if ( l < L ) {
    return Lattice::spindex_type::up;
  } else {
    return Lattice::spindex_type::down;
  }
}


vector<Lattice::spindex> Lattice::get_Xnn( Lattice::spindex l, unsigned int X ) const
{
  vector<Lattice::spindex> Xnn;
  get_Xnn( l, X, &Xnn );
  return Xnn;
}
