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

#ifndef MCCRESULTS_H_INCLUDED
#define MCCRESULTS_H_INCLUDED

#include <vector>
#include <numeric>

#include <boost/optional.hpp>

#include <eigen3/Eigen/Core>

#include "fptype.hpp"


template<typename T>
struct UncertainQuantity {

  T mean, sigma;

  UncertainQuantity() {}

  UncertainQuantity( T mean_init, T sigma_init )
    : mean( mean_init ), sigma( sigma_init ) { }

  UncertainQuantity( const std::vector<T>& binmeans ) {
    // calculate the average of the binmeans
    mean =
      accumulate( binmeans.begin(), binmeans.end(), 0.f ) /
      static_cast<fptype>( binmeans.size() );

    // calculate the variance of the binmeans
    fptype binmeans_variance =
      static_cast<fptype>( binmeans.size() ) /
      static_cast<fptype>( binmeans.size() - 1 ) *
      (
        accumulate(
          binmeans.begin(), binmeans.end(), 0.f,
          []( fptype sum, fptype m ) { return sum + m * m; }
        ) / static_cast<fptype>( binmeans.size() )
        - mean * mean
      );

    // uncertainty of the mean is sqrt(variance / num_bins)
    sigma = sqrt( binmeans_variance / static_cast<fptype>( binmeans.size() ) );
  }
};


struct MCCResults {

  bool success;

  boost::optional< UncertainQuantity<fptype> > E;
  boost::optional< Eigen::VectorXfp > Deltak;
  boost::optional< Eigen::MatrixXfp > Deltak_Deltakprime;
  boost::optional< Eigen::VectorXfp > Deltak_E;
};

#endif // MCCRESULTS_H_INCLUDED
