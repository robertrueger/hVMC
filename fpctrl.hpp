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

#ifndef FPCTRL_H_INCLUDED
#define FPCTRL_H_INCLUDED

#include <eigen3/Eigen/Core>

#include "macros.h"
#include "fptype.hpp"


struct FPDevStat final {

  // the desired deviation
  cl_fptype target;

  // the number of recalculations
  cl_uint recalcs;

  // the number of misses
  cl_uint misses;

  // the number of target hits
  cl_uint hits;

  // the number of misses by more than an order of magnitude
  cl_uint mag1_misses;

  // the number of hits better than an order of magnitude
  cl_uint mag1_hits;

  FPDevStat( cl_fptype target_init )
    : target( target_init ),
      recalcs( 0 ), misses( 0 ), hits( 0 ),
      mag1_misses( 0 ), mag1_hits( 0 ) {}

  void add( cl_fptype dev );
};


cl_fptype calc_deviation(
  const Eigen::MatrixXfp& approx, const Eigen::MatrixXfp& exact
);

#endif // FPCTRL_H_INCLUDED
