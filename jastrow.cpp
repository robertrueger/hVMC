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

using namespace std;


Jastrow::Jastrow(
  const shared_ptr<Lattice>& lat_init, const Eigen::VectorXd& v_init )
  : lat( lat_init )
{
  std::set<Lattice::irrspidxrel> isrs = lat->get_all_irrspidxrels();
  assert( !isrs.empty() );

  // find maximum j in v_0j
  const unsigned int max_j_v = *( max_element( isrs.begin(), isrs.end() ) );
  // resize internal vector so that it can hold all needed Jastrows
  isr_v.resize( max_j_v + 1 );

  // delete irreducible spindex relation corresponding to the largest distance
  isrs.erase( lat->get_maxdist_irrspidxrel() );

  // find maximum j in v_0j if the largest distance v is missing
  const unsigned int max_j_vparnum = *( max_element( isrs.begin(), isrs.end() ) );
  // resize vector that maps irreducible irrspidxrel -> variational parameter number
  isr_vparnum.resize( max_j_vparnum + 1 );
  // save the total number of variational parameters
  num_vpar = isrs.size();

  // write the variational parameters from v_init to the right elements of
  // isr_v (make sure the total number is correct first)
  assert( static_cast<unsigned int>( v_init.size() ) == isrs.size() );
  unsigned int reader = 0;
  for ( auto isrs_it = isrs.begin(); isrs_it != isrs.end(); ++isrs_it ) {
    isr_v.at(       *isrs_it ) = v_init( reader );
    isr_vparnum.at( *isrs_it ) = reader;
    ++reader;
  }
  assert( reader == v_init.size() );

  // set the parameter corresponding to the largest distance to 0
  isr_v.at( lat->get_maxdist_irrspidxrel() ) = 0.0;
}



double Jastrow::operator()( Lattice::spindex i, Lattice::spindex j ) const
{
  assert( isr_v.size() > lat->reduce_spidxrel( i, j ) );
  return isr_v[ lat->reduce_spidxrel( i, j ) ];
}



double Jastrow::onsite() const
{
  assert( isr_v.size() > 0 );
  return isr_v[0];
}



void Jastrow::set( Lattice::spindex i, Lattice::spindex j, double v_new )
{
  isr_v.at( lat->reduce_spidxrel( i, j ) ) = v_new;
}



unsigned int Jastrow::get_num_vpar() const
{
  return num_vpar;
}



unsigned int Jastrow::get_vparnum( unsigned int isr ) const
{
  assert( isr_vparnum.size() > isr );
  return isr_vparnum[ isr ];
}
