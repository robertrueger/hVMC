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
#include <cmath>

#include <boost/optional.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/optional.hpp>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>


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

    // if binmeans variance is negative this is floating point precision
    // problem (cancellation when subtracting mean^2)
    binmeans_variance = std::max( binmeans_variance, 0.0 );

    // uncertainty of the mean is sqrt(variance / num_bins)
    sigma = std::sqrt( binmeans_variance / static_cast<Ti>( binmeans.size() ) );
  }

  // make UncertainQuantity serializable
  friend class boost::serialization::access;
  template<class Archive>
  void serialize( Archive& ar, const unsigned int ) {
    ar & mean;
    ar & sigma;
  }
};


struct MCCResults {

  boost::optional< UncertainQuantity<double> > E;
  boost::optional< Eigen::VectorXd > Deltak;
  boost::optional< Eigen::MatrixXd > Deltak_Deltakprime;
  boost::optional< Eigen::VectorXd > Deltak_E;
  boost::optional< UncertainQuantity<double> > dblocc;
  boost::optional< Eigen::MatrixXd > nncorr;
  boost::optional< Eigen::MatrixXd > sscorr;
  boost::optional< std::vector<Eigen::VectorXi> > pconfs;

  void write_to_files( const boost::filesystem::path& dir ) const;
};
std::ostream& operator<<( std::ostream& out, const MCCResults& res );

#endif // MCCRESULTS_H_INCLUDED
