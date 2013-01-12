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

#include "simresults.hpp"

#include <cmath>

using namespace std;


BasicSimResults combine_results( const vector<BasicSimResults>& result_list )
{
  assert( result_list.size() > 0 );

  BasicSimResults combined_results;
  combined_results.success = true;

  for ( unsigned int i = 0; i < result_list.size(); ++i ) {
    combined_results.success &= result_list[i].success;
    combined_results.E += result_list[i].E;
    combined_results.sigma_E += result_list[i].sigma_E * result_list[i].sigma_E;
    combined_results.var_E_l += result_list[i].var_E_l;
  }

  combined_results.E /= static_cast<fptype>( result_list.size() );
  combined_results.sigma_E
    = sqrt( combined_results.sigma_E ) / static_cast<fptype>( result_list.size() );
  combined_results.var_E_l /= static_cast<fptype>( result_list.size() );

  return combined_results;
}
