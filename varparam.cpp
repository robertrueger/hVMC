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

#include "varparam.hpp"

#include <iostream>
#include <fstream>
#include <set>

#include <boost/filesystem.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "serialization_eigen.hpp"
#include "mccrun_prepare.hpp"
#include "lattice.hpp"

using namespace std;
namespace fs  = boost::filesystem;
namespace ar  = boost::archive;


Eigen::VectorXfp get_initial_varparam( const Options& opts )
{
  // determine how many variational parameters there are
  unsigned int num_vpars
    = prepare_lattice( opts )->irreducible_idxrel_list().size() - 1;

  if ( fs::exists(
         opts["calc.working-dir"].as<fs::path>() /
         "opt_vpar_final.dat" )
      ) {
    // read the variational parameters from a file
    ifstream vpar_file( (
      opts["calc.working-dir"].as<fs::path>() /
       "opt_vpar_final.dat"
    ).string() );
    ar::text_iarchive vpar_archive( vpar_file );
    Eigen::VectorXfp vpar;
    vpar_archive >> vpar;
    if ( vpar.size() != num_vpars ) {
      cerr << "ERROR: variational parameter file does not have the right number "
              "of variational parameters!" << endl;
      exit( 1 );
    }
    return vpar;
  } else {
    // set all variational parameters to zero
    return Eigen::VectorXfp::Zero( num_vpars );
  }
}
