/*
 * Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
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

#include <iostream>
#include <vector>
#include <chrono>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>

#include "macros.h"
#include "options.hpp"
#include "varparam.hpp"
#include "simresults.hpp"
#include "simrun.hpp"

using namespace std;
namespace mpi = boost::mpi;


int main( int argc, char* argv[] )
{
  // initialize mpi
  mpi::environment  mpienv( argc, argv );
  mpi::communicator mpiflock;

  if ( mpiflock.rank() == 0 ) { // only flock leader prints welcome message
    cout << endl;
    cout << "    ==========================================" << endl;
    cout << "    | hVMC - hubbard Variational Monte Carlo |" << endl;
    cout << "    ==========================================" << endl;
    cout << endl;
  }

  // read options from the command line
  const Options& opts = read_options( argc, argv, mpiflock );
  // TODO: output all options if verbose (Robert Rueger, 2012-11-02 13:31)

  if ( mpiflock.rank() == 0 ) {
    // ----- thread is leading the flock -----

    // start the stopwatch
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

    if ( mpiflock.size() != 1 ) {
      cout << ":: Broadcasting initial variational parameters" << endl;
    }
    VariationalParameters vpar = get_initial_varparam( opts );
    mpi::broadcast( mpiflock, vpar, 0 );
    if ( mpiflock.size() != 1 ) {
      cout << ":: Starting the simulation" << endl;
    }

    vector<BasicSimResults> flockresults;

    cout << endl;
    cout << "----- SIMULATION OUTPUT ------------------------------" << endl;
    const BasicSimResults& myresult = simrun_basic( opts, vpar, mpiflock );
    cout << "------------------------------------------------------" << endl;
    cout << endl;

    if ( mpiflock.size() != 1 ) {
      cout << ":: Collecting the flock's results" << endl;
    }
    mpi::gather( mpiflock, myresult, flockresults, 0 );

    // stop the stopwatch and calculate the elapsed time
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    const double total_time = chrono::duration<double>( t2 - t1 ).count();
    const double time_per_mcs = total_time / static_cast<double>
                                ( opts["sim.num-bins"].as<unsigned int>() *
                                  opts["sim.num-binmcs"].as<unsigned int>() );
    cout << ":: Simulation has finished in " << total_time << " sec" << endl;
    cout << "   -> Total performance = "
         << 1.0 / time_per_mcs << " effMCS/sec" << endl;

    // display results of the threads individually
    if ( opts.count( "verbose" ) ) {
      cout << ":: Displaying individual results" << endl;
      for ( unsigned int i = 0; i < flockresults.size(); ++i ) {
        cout << "Result of thread " << i << endl;
        cout << "       E = " << flockresults[i].E << endl;
        cout << " sigma_E = " << flockresults[i].sigma_E << endl;
        cout << "var(E_l) = " << flockresults[i].var_E_l << endl;
      }
    }

    const BasicSimResults& combined_results = combine_results( flockresults );
    if ( mpiflock.size() == 1 ) {
      cout << ":: Simulation results" << endl;
    } else {
      cout << ":: Combined result" << endl;
    }
    cout << "       E = " << combined_results.E << endl;
    cout << " sigma_E = " << combined_results.sigma_E << endl;
    cout << "var(E_l) = " << combined_results.var_E_l << endl;
    cout << endl;

  } else {
    // ----- thread is following the flock -----

    // receive variational parameters from the leader
    VariationalParameters vpar;
    mpi::broadcast( mpiflock, vpar, 0 );

    // run the simulation
    const BasicSimResults& myresult = simrun_basic( opts, vpar, mpiflock );

    // send results to the leader
    mpi::gather( mpiflock, myresult, 0 );
  }

  return 0;
}
