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

#include "sched.hpp"

#include <set>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/set.hpp>

#include "serialization_eigen.hpp"
#include "mccresults.hpp"
#include "mccrun.hpp"
#include "msgtags.hpp"
#include "fptype.hpp"
#include "observables.hpp"
#include "varparam.hpp"

using namespace std;
namespace mpi = boost::mpi;


void sched_master( const Options& opts, const mpi::communicator& mpicomm )
{
  schedmsg_t schedmsg;

  // TODO: use options to decide what to do (Robert Rueger, 2013-02-20 12:12)
  // for now: just run a single Monte Carlo cycle ...

  // tell everyone that we want to start a Monte Carlo cycle
  schedmsg = SCHEDMSG_START_MCC;
  mpi::broadcast( mpicomm, schedmsg, 0 );

  // prepare the input data for the Monte Carlo Cycle
  const Eigen::VectorXfp vpar = get_initial_varparam( opts );

  set<observables_t> obs;
  obs.insert( OBSERVABLE_E );

  // send the varparams and the set of observables to the slaves
  mpi::broadcast( mpicomm, vpar, 0);
  mpi::broadcast( mpicomm, obs,  0);

  // run master part of the Monte Carlo cycle
  const MCCResults& res
    = mccrun_master( opts, vpar, obs, mpicomm );

  // everything done, tell everyone to quit!
  schedmsg = SCHEDMSG_EXIT;
  mpi::broadcast( mpicomm, schedmsg, 0 );
}


void sched_slave( const Options& opts, const mpi::communicator& mpicomm )
{
  schedmsg_t schedmsg;
  do {
    mpi::broadcast( mpicomm, schedmsg, 0 );

    if ( schedmsg == SCHEDMSG_START_MCC ) {

      // get variational parameters and set of observables from master
      Eigen::VectorXfp vpar;
      mpi::broadcast( mpicomm, vpar, 0);
      set<observables_t> obs;
      mpi::broadcast( mpicomm, obs,  0);

      // run slave part of the Monte Carlo cycle
      mccrun_slave( opts, vpar, obs, mpicomm );

    }

  } while ( schedmsg != SCHEDMSG_EXIT );
}
