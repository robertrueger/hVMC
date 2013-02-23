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

#ifndef FPCTRL_H_INCLUDED
#define FPCTRL_H_INCLUDED

#include <boost/serialization/access.hpp>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "macros.h"
#include "fptype.hpp"


struct FPDevStat {

  // the desired deviation
  fptype target;

  // the number of recalculations
  unsigned int recalcs;

  // the number of misses
  unsigned int misses;

  // the number of target hits
  unsigned int hits;

  // the number of misses by more than an order of magnitude
  unsigned int mag1_misses;

  // the number of hits better than an order of magnitude
  unsigned int mag1_hits;

  FPDevStat()
    : target( 0.f ),
      recalcs( 0 ), misses( 0 ), hits( 0 ),
      mag1_misses( 0 ), mag1_hits( 0 ) { }

  FPDevStat( fptype target_init )
    : target( target_init ),
      recalcs( 0 ), misses( 0 ), hits( 0 ),
      mag1_misses( 0 ), mag1_hits( 0 ) { }

  void add( fptype dev );

  // make FPDevStat serializable
  friend class boost::serialization::access;
  template<class Archive>
  void serialize( Archive& ar, const unsigned int ) {
    ar & target;
    ar & recalcs;
    ar & misses;
    ar & hits;
    ar & mag1_misses;
    ar & mag1_hits;
  }
};

FPDevStat operator+( const FPDevStat& lhs, const FPDevStat& rhs );

fptype calc_deviation(
  const Eigen::MatrixXfp& approx, const Eigen::MatrixXfp& exact
);

#endif // FPCTRL_H_INCLUDED
