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

#include "obs_dblocc.hpp"

#include <numeric>

#include <boost/mpi/collectives.hpp>

using namespace std;
namespace mpi = boost::mpi;


ObservableDoubleOccupancy::ObservableDoubleOccupancy()
  : Observable( OBSERVABLE_DOUBLE_OCCUPANCY_DENSITY ) { }


void ObservableDoubleOccupancy::measure(
  const HubbardModelVMC& model, ObservableCache& cache )
{
  if ( !cache.dblocc ) {
    cache.dblocc = model.dblocc_dens();
  }

  dblocc_currentbin.push_back( cache.dblocc.get() );
}


void ObservableDoubleOccupancy::completebin()
{
  dblocc_binmeans.push_back(
    accumulate( dblocc_currentbin.begin(), dblocc_currentbin.end(), 0.0 ) /
    static_cast<double>( dblocc_currentbin.size() )
  );
  dblocc_currentbin.clear();
}


void ObservableDoubleOccupancy::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );
  vector< vector<double> > binmeans_collector;
  mpi::gather( mpicomm, dblocc_binmeans, binmeans_collector, 0 );

  vector<double> dblocc_binmeans_all;
  for ( auto it = binmeans_collector.begin();
        it != binmeans_collector.end();
        ++it ) {
    dblocc_binmeans_all.insert(
        dblocc_binmeans_all.end(), it->begin(), it->end()
    );
  }

  results.dblocc = UncertainQuantity<double>( dblocc_binmeans_all );
}


void ObservableDoubleOccupancy::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );
  mpi::gather( mpicomm, dblocc_binmeans, 0 );
}
