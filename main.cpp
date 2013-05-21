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

#include <iostream>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/filesystem.hpp>

#include "macros.h"
#include "options.hpp"
#include "sched.hpp"

using namespace std;
namespace mpi = boost::mpi;
namespace fs  = boost::filesystem;


int main( int argc, char* argv[] )
{
  // initialize mpi
  mpi::environment  mpienv( argc, argv );
  mpi::communicator mpicomm;

  if ( mpicomm.rank() == 0 ) { // only master prints welcome message
    cout << endl;
    cout << "    ==========================================" << endl;
    cout << "    | hVMC - hubbard Variational Monte Carlo |" << endl;
    cout << "    ==========================================" << endl;
    cout << endl;
  }

  // read options from the command line
  const Options& opts = read_options( argc, argv, mpicomm.rank() == 0 );
  // TODO: output all options if verbose (Robert Rueger, 2012-11-02 13:31)

  if ( mpicomm.rank() == 0 ) { // this process is the master
    // create output directory
    fs::create_directory( opts["calc.working-dir"].as<fs::path>() );

    sched_master( opts, mpicomm );
  } else { // this process is a slave
    sched_slave( opts, mpicomm );
  }

  return 0;
}
