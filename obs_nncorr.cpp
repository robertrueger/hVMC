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

#include "obs_nncorr.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <boost/mpi/collectives.hpp>

#include "macros.h"
#include "serialization_eigen.hpp"

using namespace std;
namespace mpi = boost::mpi;


ObservableDensityDensityCorrelation::ObservableDensityDensityCorrelation(
  unsigned int L )
  : Observable( OBSERVABLE_DENSITY_DENSITY_CORRELATION ),
    thisbin_nncorr_sum(
      Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic>::Zero( L, L )
    ),
    thisbin_count( 0 ),
    binmean_nncorr_sum( Eigen::MatrixXd::Zero( L, L ) ),
    binmean_count( 0 ) { }


void ObservableDensityDensityCorrelation::measure(
  const HubbardModelVMC& model, ObservableCache& cache )
{
  if ( !cache.n ) {
    cache.n = model.n();
  }
  const Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>& n_current = cache.n.get();

  thisbin_nncorr_sum += n_current * n_current.transpose();
  ++thisbin_count;

#if VERBOSE >= 1
  cout << "ObservableDensityDensityCorrelation::measure() : thisbin_nncorr_sum = "
       << endl << thisbin_nncorr_sum << endl;
#endif
}


void ObservableDensityDensityCorrelation::completebin()
{
  binmean_nncorr_sum
    += thisbin_nncorr_sum.cast<double>() / static_cast<double>( thisbin_count );
  ++binmean_count;

  thisbin_nncorr_sum.setZero();
  thisbin_count = 0;
}


void ObservableDensityDensityCorrelation::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );

  vector<Eigen::MatrixXd> binmeans_collector;
  mpi::gather( mpicomm, binmean_nncorr_sum, binmeans_collector, 0 );
  vector<unsigned int> binmeans_collector_count;
  mpi::gather( mpicomm, binmean_count, binmeans_collector_count, 0 );

  results.nncorr
    = accumulate(
        binmeans_collector.begin() + 1,
        binmeans_collector.end(),
        binmeans_collector.front()
      ) / static_cast<double>(
        accumulate(
          binmeans_collector_count.begin(),
          binmeans_collector_count.end(),
          0
        )
      );
}


void ObservableDensityDensityCorrelation::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );

  mpi::gather( mpicomm, binmean_nncorr_sum, 0 );
  mpi::gather( mpicomm, binmean_count, 0 );
}
