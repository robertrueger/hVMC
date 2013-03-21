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


ObservableDensityDensityCorrelation::ObservableDensityDensityCorrelation()
  : Observable( OBSERVABLE_DENSITY_DENSITY_CORRELATION ),
    this_bin_num_measurements( 0 ) { }


void ObservableDensityDensityCorrelation::measure(
  HubbardModelVMC& model, ObservableCache& cache )
{
  if ( !cache.n ) {
    cache.n = model.n();
  }

  const Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic>&
    n_current = cache.n.get();

  if ( nncorr_sum.size() == 0 ) {
    // first use of nncorr_sum
    nncorr_sum.setZero( n_current.size(), n_current.size() );
  } else {
    assert( nncorr_sum.cols() == n_current.size() );
    assert( nncorr_sum.rows() == n_current.size() );
  }

  nncorr_sum += n_current * n_current.transpose();
  ++this_bin_num_measurements;

#if VERBOSE >= 1
  cout << "ObservableDensityDensityCorrelation::measure() : nncorr_sum = " << endl
       << nncorr_sum << endl;
#endif
}


void ObservableDensityDensityCorrelation::completebin()
{
  nncorr_binmeans.push_back(
    nncorr_sum.cast<fptype>() / static_cast<fptype>( this_bin_num_measurements )
  );

#if VERBOSE >= 1
  cout << "ObservableDensityDensityCorrelation::completebin() : binmean = " << endl
       << *( nncorr_binmeans.rbegin() ) << endl;
#endif

  nncorr_sum.setZero();
  this_bin_num_measurements = 0;
}


void ObservableDensityDensityCorrelation::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );
  vector< vector<Eigen::MatrixXfp> > binmeans_collector;
  mpi::gather( mpicomm, nncorr_binmeans, binmeans_collector, 0 );

  vector< Eigen::MatrixXfp > nncorr_binmeans_all;
  for ( auto it = binmeans_collector.begin();
        it != binmeans_collector.end();
        ++it ) {
    nncorr_binmeans_all.insert( nncorr_binmeans_all.end(), it->begin(), it->end() );
  }
  assert( !nncorr_binmeans_all.empty() );


  Eigen::MatrixXfp nncorr_binmeans_all_sum;
  nncorr_binmeans_all_sum.setZero(
    nncorr_binmeans_all[0].rows(),
    nncorr_binmeans_all[0].cols()
  );
  for ( auto it = nncorr_binmeans_all.begin();
        it != nncorr_binmeans_all.end();
        ++it ) {
    nncorr_binmeans_all_sum += *it;
#if VERBOSE >= 2
    cout << it->transpose() << endl;
#endif
  }
  results.nncorr
    = nncorr_binmeans_all_sum / static_cast<fptype>( nncorr_binmeans_all.size() );
}


void ObservableDensityDensityCorrelation::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );
  mpi::gather( mpicomm, nncorr_binmeans, 0 );
}
