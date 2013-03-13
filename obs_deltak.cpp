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

#include <boost/mpi/collectives.hpp>

using namespace std;
namespace mpi = boost::mpi;


ObservableDeltaK::ObservableDeltaK() : Observable( OBSERVABLE_DELTAK ) { }


void ObservableDeltaK::measure( HubbardModelVMC& model )
{
  model.Delta_k();
}


void ObservableDeltaK::completebin()
{

}


void ObservableDeltaK::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );

}


void ObservableDeltaK::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );

}
