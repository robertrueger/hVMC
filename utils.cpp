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

#include "utils.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

using namespace std;


void ostream_setup( ostream& stream )
{
  // stream << setiosflags( ios::scientific );
  stream.setf( ios::showpos );
  stream.precision( numeric_limits<double>::digits10 + 1 );
}


unsigned int uintsqrt( unsigned int n ) {
 return static_cast<unsigned int>(
   std::floor( std::sqrt( static_cast<double>( n ) ) + 0.5 )
 );
}


bool is_perfect_square( unsigned int n )
{
  unsigned int t = uintsqrt( n );
  return t * t == n;
}
