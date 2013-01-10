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

#if VERBOSE >= 1
# include <iostream>
#endif

#include <set>
#include <algorithm>
#include <cmath>

using namespace std;


Jastrow::Jastrow( Lattice* lat_init, const vector<fptype>& v_init )
  : lat( lat_init )
{
  const std::set<unsigned int>& irr_idxrels = lat->irreducible_idxrel_list();
  assert( !irr_idxrels.empty() );

  // find maximum j in v_0j
  const unsigned int max_j =
    *( max_element( irr_idxrels.begin(), irr_idxrels.end() ) );

  // resize internal vector so that it can hold all needed Jastrows
  idxrel_expv.resize( max_j + 1 );

  // write the variational parameters from v_init to the right elements of
  // idxrel_expv (make sure the total number is correct first)
  assert( v_init.size() == irr_idxrels.size() );
  unsigned int reader = 0;
  for ( auto irr_idxrel_it = irr_idxrels.begin();
        irr_idxrel_it != irr_idxrels.end();
        ++irr_idxrel_it ) {
    idxrel_expv.at( *irr_idxrel_it ) = std::exp( v_init.at( reader ) );
    ++reader;
  }
  assert( reader == v_init.size() );
}



void Jastrow::randomize( fptype min, fptype max, mt19937* rng )
{
  const std::set<unsigned int>& irr_idxrels = lat->irreducible_idxrel_list();

  for ( auto it = irr_idxrels.begin(); it != irr_idxrels.end(); ++it ) {
    idxrel_expv[*it]
      = std::exp( uniform_real_distribution<fptype>( min, max )( *rng ) );
  }
}



fptype Jastrow::operator()( unsigned int i, unsigned int j ) const
{
  assert( idxrel_expv.size() > lat->reduce_idxrel( i, j ) );
  return std::log( idxrel_expv[ lat->reduce_idxrel( i, j ) ] );
}



fptype Jastrow::exp( unsigned int i, unsigned int j ) const
{
  assert( idxrel_expv.size() > lat->reduce_idxrel( i, j ) );
  return idxrel_expv[ lat->reduce_idxrel( i, j ) ];
}



fptype Jastrow::exp_onsite() const
{
  assert( idxrel_expv.size() > 0 );
  return idxrel_expv[0];
}



void Jastrow::set( unsigned int i, unsigned int j, fptype v_new )
{
  idxrel_expv.at( lat->reduce_idxrel( i, j ) ) = std::exp( v_new );
}
