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

#include "jastrow.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <set>
#include <algorithm>
#include <cmath>

using namespace std;


Jastrow::Jastrow(
  const shared_ptr<Lattice>& lat_init, const Eigen::VectorXd& v_init )
  : lat( lat_init )
{
  std::set<unsigned int> irr_idxrels = lat->irreducible_idxrel_list();
  assert( !irr_idxrels.empty() );

  // find maximum j in v_0j
  const unsigned int max_j_expv =
    *( max_element( irr_idxrels.begin(), irr_idxrels.end() ) );
  // resize internal vector so that it can hold all needed Jastrows
  idxrel_expv.resize( max_j_expv + 1 );

  // delete irreducible index relation corresponding to the largest distance
  irr_idxrels.erase( lat->irreducible_idxrel_maxdist() );

  // find maximum j in v_0j if the largest distance v is missing
  const unsigned int max_j_vparnum =
    *( max_element( irr_idxrels.begin(), irr_idxrels.end() ) );
  // resize vector: irreducible idxrel -> variational parameter number
  idxrel_vparnum.resize( max_j_vparnum + 1 );
  // save the total number of variational parameters
  num_vpar = irr_idxrels.size();

  // write the variational parameters from v_init to the right elements of
  // idxrel_expv (make sure the total number is correct first)
  assert( static_cast<unsigned int>( v_init.size() ) == irr_idxrels.size() );
  unsigned int reader = 0;
  for ( auto irr_idxrel_it = irr_idxrels.begin();
        irr_idxrel_it != irr_idxrels.end();
        ++irr_idxrel_it ) {
    idxrel_expv.at(    *irr_idxrel_it ) = std::exp( v_init( reader ) );
    idxrel_vparnum.at( *irr_idxrel_it ) = reader;
    ++reader;
  }
  assert( reader == v_init.size() );

  // set the parameter corresponding to the largest distance to 1
  idxrel_expv.at( lat->irreducible_idxrel_maxdist() ) = 1.f;
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



void Jastrow::set( unsigned int i, unsigned int j, double v_new )
{
  idxrel_expv.at( lat->reduce_idxrel( i, j ) ) = std::exp( v_new );
}



unsigned int Jastrow::get_num_vpar() const
{
  return num_vpar;
}



unsigned int Jastrow::get_vparnum( unsigned int irr_idxrel ) const
{
  assert( idxrel_vparnum.size() > irr_idxrel );
  return idxrel_vparnum[ irr_idxrel ];
}
