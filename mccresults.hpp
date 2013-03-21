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
#include <iosfwd>

#include <boost/optional.hpp>
#include <boost/filesystem/path.hpp>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "fptype.hpp"


template <typename Ti>
struct UncertainQuantity {

  Ti mean, sigma;

  UncertainQuantity() {}

  UncertainQuantity( Ti mean_init, Ti sigma_init )
    : mean( mean_init ), sigma( sigma_init ) { }

  template< typename Te >
  UncertainQuantity( const std::vector<Te>& binmeans ) {
    // calculate the average of the binmeans
    mean =
      static_cast<Ti>(
        accumulate( binmeans.begin(), binmeans.end(), Te(0) )
      ) / static_cast<Ti>( binmeans.size() );

    // calculate the variance of the binmeans
    Ti binmeans_variance =
      static_cast<Ti>( binmeans.size() ) /
      static_cast<Ti>( binmeans.size() - 1 ) *
      (
        static_cast<Ti>( accumulate(
          binmeans.begin(), binmeans.end(), Te(0),
          []( Te sum, Te m ) { return sum + m * m; }
        ) ) / static_cast<Ti>( binmeans.size() )
        - mean * mean
      );

    // uncertainty of the mean is sqrt(variance / num_bins)
    sigma = sqrt( binmeans_variance / static_cast<Ti>( binmeans.size() ) );
  }
};


struct MCCResults {

  bool success;

  boost::optional< UncertainQuantity<fptype> > E;
  boost::optional< Eigen::VectorXfp > Deltak;
  boost::optional< Eigen::MatrixXfp > Deltak_Deltakprime;
  boost::optional< Eigen::VectorXfp > Deltak_E;
  boost::optional< UncertainQuantity<fptype> > dblocc;
  boost::optional< Eigen::MatrixXfp > nncorr;

  void write_to_files( const boost::filesystem::path& dir ) const;
};
std::ostream& operator<<( std::ostream& out, const MCCResults& res );

#endif // MCCRESULTS_H_INCLUDED
