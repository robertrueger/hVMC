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

#include <boost/program_options.hpp>

#include "macros.h"
#include "options.hpp"
#include "varparam.hpp"
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

  // read options from the command line
  const Options& opts = read_options( argc, argv );
  // TODO: output all options if verbose (Robert Rueger, 2012-11-02 13:31)

  // initialize the variational parameters
  const VariationalParameters& vpar = get_initial_varparam( opts );

  // run the simulation
  const BasicSimResults& results = simrun_basic( opts, vpar );

  // output results
  cout << ":: Results" << endl;
  cout << "       E = " << results.E << endl;
  cout << " sigma_E = " << results.sigma_E << endl;
  cout << "var(E_l) = " << results.var_E_l << endl;

  return 0;
}
