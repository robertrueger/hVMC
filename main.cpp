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

#include <eigen3/Eigen/Core>

#include <boost/program_options.hpp>

#include "macros.h"
#include "utils.hpp"
#include "options.hpp"
#include "simresults.hpp"
#include "simrun.hpp"

using namespace std;
namespace po = boost::program_options;


int main( int argc, char const* argv[] )
{
  cout << endl;
  cout << "    ==========================================" << endl;
  cout << "    | hVMC - hubbard Variational Monte Carlo |" << endl;
  cout << "    ==========================================" << endl;
  cout << endl;

  // make cout and cerr output numbers in scientific notation
  // ostream_setup( cout );
  // ostream_setup( cerr );

  // initialize Eigen for parallel execution
  // (not sure if that is necessary, but better safe than sorry)
  Eigen::initParallel(); 

  // read options from the command line
  const Options& opts = read_options( argc, argv );

  // TODO: output all options if verbose (Robert Rueger, 2012-11-02 13:31)
 
  // run the simulation
  const BasicSimResults& results = simrun_basic( opts );

  // output results
  cout << ":: Results" << endl;
  cout << "       E = " << results.E << endl;
  cout << " sigma_E = " << results.sigma_E << endl;
  cout << "var(E_l) = " << results.var_E_l << endl;

  return 0;
}
