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

#include "obs_pconfs.hpp"

#include <boost/mpi/collectives.hpp>

#include "macros.h"
#include "serialization_eigen.hpp"

using namespace std;
namespace mpi = boost::mpi;


ObservableParticleConfigurations::ObservableParticleConfigurations()
  : thisbin_recorded( false ) { }


void ObservableParticleConfigurations::measure(
  const HubbardModelVMC& model, ObservableCache& )
{
  if ( !thisbin_recorded ) {
    site_occs.push_back( model.particleconf() );
    thisbin_recorded = true;
  }
}


void ObservableParticleConfigurations::completebin()
{
  thisbin_recorded = false;
}


void ObservableParticleConfigurations::collect_and_write_results(
  const boost::mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );
  vector< vector<Eigen::VectorXi> > site_occs_collector;
  mpi::gather( mpicomm, site_occs, site_occs_collector, 0 );

  vector<Eigen::VectorXi> site_occs_all;
  for ( auto it = site_occs_collector.begin();
        it != site_occs_collector.end();
        ++it ) {
    site_occs_all.insert( site_occs_all.end(), it->begin(), it->end() );
  }

  results.pconfs = site_occs_all;
}


void ObservableParticleConfigurations::send_results_to_master(
  const boost::mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );
  mpi::gather( mpicomm, site_occs, 0 );
}
