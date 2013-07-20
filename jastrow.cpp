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
  std::set<Lattice::irridxrel> iirs = lat->get_all_irridxrels();
  assert( !iirs.empty() );

  // find maximum j in v_0j
  const unsigned int max_j_v = *( max_element( iirs.begin(), iirs.end() ) );
  // resize internal vector so that it can hold all needed Jastrows
  iir_v.resize( max_j_v + 1 );

  // delete irreducible index relation corresponding to the largest distance
  iirs.erase( lat->get_maxdist_irridxrel() );

  // find maximum j in v_0j if the largest distance v is missing
  const unsigned int max_j_vparnum = *( max_element( iirs.begin(), iirs.end() ) );
  // resize vector that maps irreducible irridxrel -> variational parameter number
  iir_vparnum.resize( max_j_vparnum + 1 );
  // save the total number of variational parameters
  num_vpar = iirs.size();

  // write the variational parameters from v_init to the right elements of
  // iir_v (make sure the total number is correct first)
  assert( static_cast<unsigned int>( v_init.size() ) == iirs.size() );
  unsigned int reader = 0;
  for ( auto iirs_it = iirs.begin(); iirs_it != iirs.end(); ++iirs_it ) {
    iir_v.at(       *iirs_it ) = v_init( reader );
    iir_vparnum.at( *iirs_it ) = reader;
    ++reader;
  }
  assert( reader == v_init.size() );

  // set the parameter corresponding to the largest distance to 0
  iir_v.at( lat->get_maxdist_irridxrel() ) = 0.0;
}



double Jastrow::operator()( Lattice::spindex i, Lattice::spindex j ) const
{
  assert( iir_v.size() > lat->reduce_idxrel( i, j ) );
  return iir_v[ lat->reduce_idxrel( i, j ) ];
}



double Jastrow::onsite() const
{
  assert( iir_v.size() > 0 );
  return iir_v[0];
}



void Jastrow::set( Lattice::spindex i, Lattice::spindex j, double v_new )
{
  iir_v.at( lat->reduce_idxrel( i, j ) ) = v_new;
}



unsigned int Jastrow::get_num_vpar() const
{
  return num_vpar;
}



unsigned int Jastrow::get_vparnum( unsigned int iir ) const
{
  assert( iir_vparnum.size() > iir );
  return iir_vparnum[ iir ];
}
