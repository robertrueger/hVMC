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

#include "varparam.hpp"

#include <set>

#include "lattice.hpp"
#include "lattice_1dchain.hpp"
#include "lattice_2dsquare.hpp"

using namespace std;


VariationalParameters get_initial_varparam( const Options& opts )
{
  VariationalParameters par;

  // set t, t' and t" as variational parameters of the determinantal part
  // TODO: generalize to other detWF (Robert Rueger, 2013-01-09 18:36)
  par.determinantal.resize(3);
  par.determinantal[0] = opts["phys.nn-hopping"].as<fptype>();
  par.determinantal[1] = opts["phys.2nd-nn-hopping"].as<fptype>();
  par.determinantal[2] = opts["phys.3rd-nn-hopping"].as<fptype>();

  // find out how many Jastrow parameters we need for this lattice
  // ... and set all of them to zero
  Lattice* lat;
  if ( opts["phys.lattice"].as<lattice_t>() == LATTICE_1DCHAIN ) {
    lat = new Lattice1DChain(
      opts["phys.num-lattice-sites"].as<unsigned int>() );
  } else {
    lat = new Lattice2DSquare(
      opts["phys.num-lattice-sites"].as<unsigned int>() );
  }
  par.jastrow.resize( lat->irreducible_idxrel_list().size(), 0.f );
  delete lat;
  lat = nullptr;

  return par;
}
