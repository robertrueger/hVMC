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

#include "jastrow.hpp"
using namespace std;


Jastrow::Jastrow( Lattice* lat_init )
  : lat( lat_init )
{
  const std::set<IrreducibleIdxRel>& irr_idxrels
    = lat->irreducible_idxrel_list();
  for ( auto it = irr_idxrels.begin(); it != irr_idxrels.end(); ++it ) {
    idxrel_v_map.insert( pair< IrreducibleIdxRel, fptype>( *it, 0.f ) );
  }
}



void Jastrow::randomize( fptype min, fptype max, mt19937* rng )
{
  for (auto it = idxrel_v_map.begin(); it != idxrel_v_map.end(); ++it ) {
    (*it).second = uniform_real_distribution<fptype>(min, max)(*rng);
  }
}



fptype Jastrow::operator()( unsigned int i, unsigned int j ) const
{
  const auto& it = idxrel_v_map.find( lat->reduce_idxrel( i, j ) );
  assert( it != idxrel_v_map.end() );
  return (*it).second;
}



void Jastrow::set( unsigned int i, unsigned int j, fptype v_new )
{
  const auto& it = idxrel_v_map.find( lat->reduce_idxrel( i, j ) );
  assert( it != idxrel_v_map.end() );
  (*it).second = v_new;
}
