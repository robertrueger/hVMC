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

#include "obs_dkdkp.hpp"

#include <boost/mpi/collectives.hpp>

using namespace std;
namespace mpi = boost::mpi;


ObservableDeltaKDeltaKPrime::ObservableDeltaKDeltaKPrime()
  : Observable( OBSERVABLE_DELTAK_DELTAKPRIME ) { }


void ObservableDeltaKDeltaKPrime::measure( HubbardModelVMC& model )
{

}


void ObservableDeltaKDeltaKPrime::completebin()
{

}


void ObservableDeltaKDeltaKPrime::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );

}


void ObservableDeltaKDeltaKPrime::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );

}
