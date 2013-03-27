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

#include "obs_deltak.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <boost/mpi/collectives.hpp>

#include "macros.h"
#include "serialization_eigen.hpp"

using namespace std;
namespace mpi = boost::mpi;


ObservableDeltaK::ObservableDeltaK()
  : Observable( OBSERVABLE_DELTAK ),
    this_bin_num_measurements( 0 ) { }


void ObservableDeltaK::measure(
  const HubbardModelVMC& model, ObservableCache& cache )
{
  if ( !cache.DeltaK ) {
    cache.DeltaK = model.Delta_k();
  }

  const Eigen::VectorXd& Dk_current = cache.DeltaK.get();

  if ( Dk_sum.size() == 0 ) {
    // first use of Dk_sum
    Dk_sum.setZero( Dk_current.size() );
  } else {
    assert( Dk_sum.size() == Dk_current.size() );
  }

  Dk_sum += Dk_current;
  ++this_bin_num_measurements;

#if VERBOSE >= 1
  cout << "ObservableDeltaK::measure() : Dk_sum = " << endl
       << Dk_sum.transpose() << endl;
#endif
}


void ObservableDeltaK::completebin()
{
  Dk_binmeans.push_back(
    Dk_sum / static_cast<double>( this_bin_num_measurements )
  );

#if VERBOSE >= 1
  cout << "ObservableDeltaK::completebin() : binmean = " << endl
       << Dk_binmeans.rbegin()->transpose() << endl;
#endif

  Dk_sum.setZero();
  this_bin_num_measurements = 0;
}


void ObservableDeltaK::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );
  vector< vector<Eigen::VectorXd> > binmeans_collector;
  mpi::gather( mpicomm, Dk_binmeans, binmeans_collector, 0 );

  vector< Eigen::VectorXd > Dk_binmeans_all;
  for ( auto it = binmeans_collector.begin();
        it != binmeans_collector.end();
        ++it ) {
    Dk_binmeans_all.insert( Dk_binmeans_all.end(), it->begin(), it->end() );
  }
  assert( !Dk_binmeans_all.empty() );


  Eigen::VectorXd Dk_binmeans_all_sum;
  Dk_binmeans_all_sum.setZero( Dk_binmeans_all[0].size() );
  for ( auto it = Dk_binmeans_all.begin(); it != Dk_binmeans_all.end(); ++it ) {
    Dk_binmeans_all_sum += *it;
#if VERBOSE >= 2
    cout << it->transpose() << endl;
#endif
  }
  results.Deltak
    = Dk_binmeans_all_sum / static_cast<double>( Dk_binmeans_all.size() );
}


void ObservableDeltaK::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );
  mpi::gather( mpicomm, Dk_binmeans, 0 );
}
