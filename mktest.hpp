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

#ifndef MANN_KENDALL_TEST_H_INCLUDED
#define MANN_KENDALL_TEST_H_INCLUDED

#include <vector>
#include <cmath>

#include <boost/math/special_functions/sign.hpp>


template <typename T>
double mktest( const std::vector<T>& data )
{
  int S = 0;
  const unsigned int n = data.size();
  for ( unsigned int i = 0; i < n - 1; ++i ) {
    for ( unsigned int j = i + 1; j < n; ++j ) {
      S += boost::math::sign( data[j] - data[i] );
    }
  }

  const double relS = std::abs( S / ( static_cast<double>( n * ( n - 1 ) / 2 ) ) );

  return relS;
}

#endif // MANN_KENDALL_TEST_H_INCLUDED
