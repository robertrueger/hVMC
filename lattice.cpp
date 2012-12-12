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

#include "lattice.hpp"
using namespace std;


cl_uint Lattice::get_spinup_site( cl_uint l ) const
{
  assert( l < 2 * L );
  return l >= L ? l - L : l;
}


cl_uint Lattice::get_spinlinked_site( cl_uint l ) const
{
  assert( l < 2 * L );

  if ( l < L ) {
    return l + L;
  } else {
    return l - L;
  }
}
