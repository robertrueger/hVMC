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

#include "obs_energy.hpp"

#include <numeric>

#include <boost/mpi/collectives.hpp>

using namespace std;
namespace mpi = boost::mpi;


ObservableEnergy::ObservableEnergy() : Observable( OBSERVABLE_E ) { }


void ObservableEnergy::measure(
  const HubbardModelVMC& model, ObservableCache& cache )
{
  if ( !cache.E ) {
    cache.E = model.E_l();
  }

  E_l_currentbin.push_back( cache.E.get() );
}


void ObservableEnergy::completebin()
{
  E_l_binmeans.push_back(
    accumulate( E_l_currentbin.begin(), E_l_currentbin.end(), 0.f ) /
    static_cast<fptype>( E_l_currentbin.size() )
  );
  E_l_currentbin.clear();
}


void ObservableEnergy::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );
  vector< vector<fptype> > binmeans_collector;
  mpi::gather( mpicomm, E_l_binmeans, binmeans_collector, 0 );

  vector< fptype > E_l_binmeans_all;
  for ( auto it = binmeans_collector.begin();
        it != binmeans_collector.end();
        ++it ) {
    E_l_binmeans_all.insert( E_l_binmeans_all.end(), it->begin(), it->end() );
  }

  results.E = UncertainQuantity<fptype>( E_l_binmeans_all );
}


void ObservableEnergy::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );
  mpi::gather( mpicomm, E_l_binmeans, 0 );
}
