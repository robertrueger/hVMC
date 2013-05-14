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
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "serialization_eigen.hpp"
#include "mccrun_prepare.hpp"
#include "lattice.hpp"
#include "detwf.hpp"

using namespace std;
namespace fs  = boost::filesystem;
namespace ar  = boost::archive;


unsigned int get_num_vpars( const Options& opts )
{
  return 7 + prepare_lattice( opts )->get_all_irridxrels().size() - 1;
}


Eigen::VectorXd get_initial_varparam( const Options& opts )
{
  // determine how many variational parameters there are
  unsigned int num_vpars = get_num_vpars( opts );

  // prepare vector
  Eigen::VectorXd init_vpar( num_vpars );

  if ( fs::exists( opts["phys.vpar-file"].as<fs::path>() ) ) {

    // read the variational parameters from a file
    ifstream vpar_file( ( opts["phys.vpar-file"].as<fs::path>() ).string() );
    ar::text_iarchive vpar_archive( vpar_file );
    Eigen::VectorXd vpar_fromfile;
    vpar_archive >> vpar_fromfile;
    if ( vpar_fromfile.size() != num_vpars ) {
      cerr << "ERROR: variational parameter file does not have the right number "
           "of variational parameters!" << endl;
      exit( 1 );
    }
    init_vpar = vpar_fromfile;

  } else {

    // choose a reasonable default value
    init_vpar.setZero();

    // set the variational t to the t in the Hubbard Hamiltonian
    init_vpar( 0 ) = opts["phys.2nd-nn-hopping"].as<double>();
    init_vpar( 1 ) = opts["phys.3rd-nn-hopping"].as<double>();

  }

  // apply variational parameter overwrites for t

  if ( opts.count( "phys.vpar-ovwrt-t2" ) ) {
    init_vpar( 0 ) = opts["phys.vpar-ovwrt-t2"].as<double>();
  }

  if ( opts.count( "phys.vpar-ovwrt-t3" ) ) {
    init_vpar( 1 ) = opts["phys.vpar-ovwrt-t3"].as<double>();
  }

  if ( !fs::exists( opts["phys.vpar-file"].as<fs::path>() ) ) {
    // set the chemical potential to a reasonable default value
    // (we have to do it here, because the ts could have been overwritten)
    const vector<double> t_vpar
      = { opts["phys.nn-hopping"].as<double>(), init_vpar( 0 ), init_vpar( 1 ) };
    init_vpar( 6 ) =
      calc_tbdetwf_chempot(
        prepare_lattice( opts ),
        opts["phys.num-electrons"].as<unsigned int>(), t_vpar
      );
  }

 // apply the other variational parameter overwrites

  if ( opts.count( "phys.vpar-ovwrt-D0" ) ) {
    init_vpar( 2 ) = opts["phys.vpar-ovwrt-D0"].as<double>();
  }

  if ( opts.count( "phys.vpar-ovwrt-D1" ) ) {
    init_vpar( 3 ) = opts["phys.vpar-ovwrt-D1"].as<double>();
  }

  if ( opts.count( "phys.vpar-ovwrt-D2" ) ) {
    init_vpar( 4 ) = opts["phys.vpar-ovwrt-D2"].as<double>();
  }

  if ( opts.count( "phys.vpar-ovwrt-D3" ) ) {
    init_vpar( 5 ) = opts["phys.vpar-ovwrt-D3"].as<double>();
  }

  if ( opts.count( "phys.vpar-ovwrt-mu" ) ) {
    init_vpar( 6 ) = opts["phys.vpar-ovwrt-mu"].as<double>();
  }

  if ( opts.count( "phys.vpar-ovwrt-J0" ) ) {
    init_vpar( 7 ) = opts["phys.vpar-ovwrt-J0"].as<double>();
  }

  return init_vpar;
}
