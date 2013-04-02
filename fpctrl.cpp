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

#include "fpctrl.hpp"


void FPDevStat::add( double dev )
{
  ++recalcs;
  if ( dev < target ) {
    ++hits;
  } else {
    ++misses;
  }
  if ( dev < 0.1 * target ) {
    ++mag1_hits;
  } else if ( dev > 10.0 * target ) {
    ++mag1_misses;
  }
}


FPDevStat operator+( const FPDevStat& lhs, const FPDevStat& rhs )
{
  assert( lhs.target == rhs.target );
  FPDevStat result( lhs.target );
  result.recalcs     = lhs.recalcs + rhs.recalcs;
  result.misses      = lhs.misses + rhs.misses;
  result.hits        = lhs.hits + rhs.hits;
  result.mag1_misses = lhs.mag1_misses + rhs.mag1_misses;
  result.mag1_hits   = lhs.mag1_hits + rhs.mag1_hits;
  return result;
}
