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
#include <chrono>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/set.hpp>

#include "serialization_eigen.hpp"
#include "mccresults.hpp"
#include "mccrun.hpp"
#include "msgtags.hpp"
#include "fptype.hpp"
#include "obs.hpp"
#include "varparam.hpp"

using namespace std;
namespace mpi = boost::mpi;


void sched_master( const Options& opts, const mpi::communicator& mpicomm )
{
  schedmsg_t schedmsg;

  // prepare the initial variational parameters
  Eigen::VectorXfp vpar = get_initial_varparam( opts );
  vpar( 0 ) = -1.f;

  // add the observables you want to measure to the set
  set<observables_t> obs;
  obs.insert( OBSERVABLE_E );
  obs.insert( OBSERVABLE_DELTAK );
  obs.insert( OBSERVABLE_DELTAK_DELTAKPRIME );
  obs.insert( OBSERVABLE_DELTAK_E );

  for ( unsigned int srit = 0; ; ++srit ) {

    cout << endl;
    cout << "====> STOCHASTIC RECONFIGURATION ITERATION " << srit << endl;
    cout << endl;

    // start the stopwatch
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

    // tell everyone that we want to start a Monte Carlo cycle
    schedmsg = SCHEDMSG_START_MCC;
    mpi::broadcast( mpicomm, schedmsg, 0 );

    // send the varparams and the set of observables to the slaves
    mpi::broadcast( mpicomm, vpar, 0 );
    mpi::broadcast( mpicomm, obs,  0 );

    // run master part of the Monte Carlo cycle
    const MCCResults& res = mccrun_master( opts, vpar, obs, mpicomm );

    // stop the stopwatch and calculate the elapsed time
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    const double total_time = chrono::duration<double>( t2 - t1 ).count();
    const double time_per_mcs = total_time / static_cast<double>
                                ( opts["sim.num-bins"].as<unsigned int>() *
                                  opts["sim.num-binmcs"].as<unsigned int>() );
    cout << ":: Simulation has finished in " << total_time << " sec" << endl;
    cout << "   Total performance = "
         << 1.0 / time_per_mcs << " effMCS/sec" << endl << endl;

    cout << ":: Simulation results" << endl << endl;
    cout << "       E = " << res.E->mean << endl;
    cout << " sigma_E = " << res.E->sigma << endl << endl;

    if ( opts.count("verbose") ) {
      cout << "Delta_k = " << endl << res.Deltak->transpose() << endl << endl;
      cout << "DkDkp = " << endl << res.Deltak_Deltakprime.get() << endl << endl;
      cout << "Dk_E = " << endl << res.Deltak_E->transpose() << endl << endl;
    }

    // calculate SR matrix and forces
    const Eigen::MatrixXfp S =
      res.Deltak_Deltakprime.get() - res.Deltak.get() * res.Deltak->transpose();
    const Eigen::VectorXfp f =
      res.Deltak.get() * res.E->mean - res.Deltak_E.get();
    const Eigen::VectorXfp dvpar =
      ( S + 0.01 * Eigen::MatrixXfp::Identity( S.rows(), S.cols() ) )
      .fullPivLu().solve( f );

    // update variational parameters
    vpar += 0.1 * dvpar;

    if ( opts.count("verbose") ) {
      cout << ":: Stochastic reconfiguration objects" << endl << endl;
      cout << "S = " << endl << S << endl << endl;
      cout << "f = " << endl << f.transpose() << endl << endl;
      cout << "dvpar = " << endl << dvpar.transpose() << endl << endl;
      cout << "vpar' = " << endl << vpar.transpose() << endl << endl;
    }
  }

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
      mpi::broadcast( mpicomm, vpar, 0 );
      set<observables_t> obs;
      mpi::broadcast( mpicomm, obs,  0 );

      // run slave part of the Monte Carlo cycle
      mccrun_slave( opts, vpar, obs, mpicomm );

    }

  } while ( schedmsg != SCHEDMSG_EXIT );
}
