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


Jastrow::Jastrow( Lattice* lat_init ) : lat( lat_init )
{
  const std::set<IrreducibleIdxRel>& irr_idxrels
    = lat->irreducible_idxrel_list();
  assert( !irr_idxrels.empty() );

  // find maximum i in v_ij (is zero on any Bravais lattice)
  const cl_uint max_i =
    max_element(
      irr_idxrels.begin(), irr_idxrels.end(),
  []( const IrreducibleIdxRel & a, const IrreducibleIdxRel & b ) {
    return a.first < b.first;
  }
    )->first;
  idxrel_expv.resize( max_i + 1 );

  // find maximum j for every possible value of i
  for ( cl_uint i = 0; i <= max_i; ++i ) {
    bool i_exists = false;
    cl_uint this_i_max_j = 0;
    for ( auto it = irr_idxrels.begin(); it != irr_idxrels.end(); ++it ) {
      if ( it->first == i ) {
        i_exists = true;
      }
      if ( it->second > this_i_max_j ) {
        this_i_max_j = it->second;
      }
    }
    if ( i_exists ) {
      idxrel_expv[i].resize( this_i_max_j + 1, 1.f );
    }
  }

#if VERBOSE >= 1
  cout << "Jastrow::Jastrow() : created the following vectors" << endl;
  for ( cl_uint i = 0; i < idxrel_expv.size(); ++i ) {
    cout << "v( " << i << ", [";
    for ( cl_uint j = 0; j < idxrel_expv[i].size(); ++j ) {
      cout << j << ",";
    }
    cout << "\b] )" << endl;
  }
#endif
}



void Jastrow::randomize( cl_fptype min, cl_fptype max, mt19937* rng )
{
  const std::set<IrreducibleIdxRel>& irr_idxrels
    = lat->irreducible_idxrel_list();

  for ( auto it = irr_idxrels.begin(); it != irr_idxrels.end(); ++it ) {
    idxrel_expv[it->first][it->second]
      = std::exp( uniform_real_distribution<cl_fptype>( min, max )( *rng ) );
  }
}



cl_fptype Jastrow::operator()( cl_uint i, cl_uint j ) const
{
  const IrreducibleIdxRel& redidx = lat->reduce_idxrel( i, j );

  assert( idxrel_expv.size() > redidx.first );
  assert( idxrel_expv[redidx.first].size() > redidx.second );

  return log( idxrel_expv[redidx.first][redidx.second] );
}



cl_fptype Jastrow::exp( cl_uint i, cl_uint j ) const
{
  const IrreducibleIdxRel& redidx = lat->reduce_idxrel( i, j );

  assert( idxrel_expv.size() > redidx.first );
  assert( idxrel_expv[redidx.first].size() > redidx.second );

  return idxrel_expv[redidx.first][redidx.second];
}



cl_fptype Jastrow::exp_onsite() const
{
  assert( idxrel_expv.size() > 0 );
  assert( idxrel_expv[0].size() > 0 );

  return idxrel_expv[0][0];
}



void Jastrow::set( cl_uint i, cl_uint j, cl_fptype v_new )
{
  const IrreducibleIdxRel& redidx = lat->reduce_idxrel( i, j );

  idxrel_expv.at( redidx.first ).at( redidx.second ) = std::exp( v_new );
}



vector<cl_fptype> Jastrow::get_reduced_raw_jastrow() const
{
  if ( idxrel_expv.size() == 1 ) {
    return idxrel_expv[0];
  } else {
    throw logic_error("the Jastrow factor is not reducible on this lattice");
  }
}
