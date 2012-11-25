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

#include "fpctrl.hpp"


void FPDevStat::add( fptype dev )
{
  ++recalcs;
  if ( dev < target ) {
    ++hits;
  } else {
    ++misses;
  }
  if ( dev < .1f * target ) {
    ++mag1_hits;
  } else if ( dev > 10.f * target ) {
    ++mag1_misses;
  }
}


fptype calc_deviation(
  const Eigen::MatrixXfp& approx, const Eigen::MatrixXfp& exact )
{
  assert( approx.size() == exact.size() );
  fptype exact_square_sum = exact.array().square().sum();
  fptype  diff_square_sum = ( approx - exact ).array().square().sum();
  return sqrt( diff_square_sum / exact_square_sum );
}
