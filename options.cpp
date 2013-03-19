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

#include "options.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>

#include "fptype.hpp"
#include "lattice.hpp"
#include "utils.hpp"

using namespace std;
namespace po  = boost::program_options;
namespace fs  = boost::filesystem;


Options read_options( int argc, char* argv[], bool is_master )
{
  Options vm;

  // define command line only options

  po::options_description clionly( "command line options" );
  clionly.add_options()
  ( "help,h", "print this help message and exit" )
  ( "version,V", "print hVMC's version and exit" )
  ( "verbose,v", "makes hVMC write additional information to stdout" )
  ( "job-file,J", po::value<fs::path>(), "job file to execute" )
  ( "output-dir,o", po::value<fs::path>()->default_value( "." ), "output directory" );
  po::positional_options_description p;
  p.add( "job-file", -1 );


  // define command line and jobfile options

  po::options_description physparam( "physical parameters" );
  physparam.add_options()

  ( "phys.nn-hopping,1",
    po::value<fptype>()->required(),
    "nearest neighbor hopping matrix element t" )

  ( "phys.2nd-nn-hopping,2",
    po::value<fptype>()->default_value( 0.f ),
    "2nd nearest neighbor hopping matrix element t'" )

  ( "phys.3rd-nn-hopping,3",
    po::value<fptype>()->default_value( 0.f ),
    "3rd nearest neighbor hopping matrix element t''" )

  ( "phys.onsite-energy,U",
    po::value<fptype>()->required(),
    "on-site energy U" )

  ( "phys.lattice,l",
    po::value<lattice_t>()->required(),
    "lattice type (1dchain, 2dsquare)" )

  ( "phys.num-lattice-sites,L",
    po::value<unsigned int>()->required(),
    "number of lattice sites" )

  ( "phys.num-electrons,N",
    po::value<unsigned int>()->required(),
    "total number of electrons" )

  ( "phys.vpars,P", po::value<fs::path>(), "variational parameter file" );

  po::options_description simset( "simulation settings" );
  simset.add_options()

  ( "sim.update-hop-maxdistance,H",
    po::value<unsigned int>()->default_value( 1 ),
    "maximum hopping distance for electronic configuration updates" )

  ( "sim.num-mcs-equil,E",
    po::value<unsigned int>()->required(),
    "number of Monte Carlo steps for equilibration" )

  ( "sim.num-bins,B",
    po::value<unsigned int>()->required(),
    "number of measurement bins" )

  ( "sim.num-binmcs,M",
    po::value<unsigned int>()->required(),
    "number of Monte Carlo steps per bin" )

  ( "sim.rng-seed,S",
    po::value<unsigned int>(),
    "random number generator seed" )

  ( "sim.sr-dt,d",
    po::value<fptype>()->default_value( 1.f ),
    "controls the SR convergence: vpar += dt * dvpar" )

  ( "sim.sr-max-refinements,R",
    po::value<unsigned int>()->default_value( 4 ),
    "number of refinements during the SR cycle" )

  ( "sim.sr-averaging-cycles,A",
    po::value<unsigned int>()->default_value( 10 ),
    "number of SR cycles to average the converged variational parameters" );

  po::options_description fpctrl( "floating point precision control" );
  fpctrl.add_options()

  ( "fpctrl.W-deviation-target",
    po::value<fptype>()->default_value( 0.001f, "0.001" ),
    "deviation target for the matrix W" )

  ( "fpctrl.W-updates-until-recalc",
    po::value<unsigned int>()
#ifdef USE_FP_DBLPREC
      ->default_value( 5000 ),
#else
      ->default_value( 500 ),
#endif
    "number of quick updates until recalculation of the matrix W" )

  ( "fpctrl.T-deviation-target",
    po::value<fptype>()->default_value( 0.001f, "0.001" ),
    "deviation target for the vector T" )

  ( "fpctrl.T-updates-until-recalc",
    po::value<unsigned int>()
#ifdef USE_FP_DBLPREC
      ->default_value( 500000 ),
#else
      ->default_value( 50000 ),
#endif
    "number of quick updates until recalculation of the vector T" );

  // define option groups for cli and jobfile
  po::options_description cmdline_options;
  cmdline_options.add( clionly ).add( physparam ).add( simset ).add( fpctrl );
  po::options_description jobfile_options;
  jobfile_options.add( physparam ).add( simset ).add( fpctrl );



  try {
    // parse options from the command line
    if ( is_master ) {
      cout << ":: Parsing command line ..." << endl;
    }
    po::store( po::command_line_parser( argc, argv ).
               options( cmdline_options ).positional( p ).run(), vm );
  } catch ( const po::error& e ) {
    if ( is_master ) {
      cerr << "Error while parsing the command line: " << e.what() << endl;
    }
    exit( 1 );
  }


  // display help or hVMC version information

  if ( vm.count( "help" ) ) {
    if ( is_master ) {
      cout << endl;
      cout << "usage: hVMC [OPTIONS] JOBFILE -o OUTDIR" << endl;
      cout << cmdline_options << endl;
    }
    exit( 0 );
  }

  if ( vm.count( "version" ) ) {
    if ( is_master ) {
      cout << endl;

      cout << "hVMC - built from git commit " << GIT_HASH << endl;

#ifndef NDEBUG
      cout << "=========== DEBUG BUILD ============" << endl;
#else
      cout << "========== RELEASE BUILD ===========" << endl;
#endif

      cout << "compiled " << __DATE__ " with ";
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#  if defined(__GNUC_PATCHLEVEL__)
      cout << "GCC " << __GNUC__ << "." << __GNUC_MINOR__
           << "." << __GNUC_PATCHLEVEL__;
#  else
      cout << "GCC " << __GNUC__ << "." << __GNUC_MINOR__;
#  endif
#elif defined(__INTEL_COMPILER)
      cout << "ICC " << __INTEL_COMPILER;
#elif defined(_MSC_VER)
      cout << "MSC " << _MSC_VER;
#else
      cout << "unknown compiler";
#endif
      cout << endl;

#ifdef USE_FP_DBLPREC
      cout << "floating point precision: double" << endl;
#else
      cout << "floating point precision: single" << endl;
#endif
#ifdef USE_CBLAS
      cout << "external CBLAS: enabled" << endl;
#else
      cout << "external CBLAS: disabled" << endl;
#endif
#ifdef EIGEN_DEFAULT_TO_ROW_MAJOR
      cout << "matrix storage order: row major" << endl;
#else
      cout << "matrix storage order: column major" << endl;
#endif
      cout << endl;

      cout
          << "Copyright (C) 2013, Robert Rueger <rueger@itp.uni-frankfurt.de>" << endl
          << "License GPLv3+: GNU GPL version 3 or later"
          " <http://gnu.org/licenses/gpl.html>" << endl
          << "This is free software: you are free to change and redistribute it." << endl
          << "There is NO WARRANTY, to the extent permitted by law." << endl << endl;
    }
    exit( 0 );
  }


  if ( vm.count( "job-file" ) ) {
    // parse the jobfile
    if ( is_master ) {
      cout << ":: Parsing jobfile ..." << endl;
    }
    ifstream jobifs( vm["job-file"].as<fs::path>().string() );
    if ( jobifs.is_open() ) {
      try {
        po::store( po::parse_config_file( jobifs, jobfile_options ), vm );
      } catch ( const po::error& e ) {
        if ( is_master ) {
          cerr << "Error while parsing the job file: " << e.what() << endl;
        }
        exit( 1 );
      }
    } else {
      if ( is_master ) {
        cerr << "Error: unable to open jobfile " << vm["job-file"].as<string>()
             << endl;
      }
      exit( 1 );
    }
  }

  try {
    // finalize the variable map
    // (will throw exception on missing options)
    po::notify( vm );
  } catch ( const po::error& e ) {
    if ( is_master ) {
      cerr << "Error in program options: " << e.what() << endl;
    }
    exit( 1 );
  }


  // check for logical errors in the physical parameters
  try {

    if ( vm["phys.num-electrons"].as<unsigned int>() % 2 != 0 ) {
      throw logic_error( "electron number must be even" );
    }

    if ( vm["phys.num-electrons"].as<unsigned int>() >
         2 * vm["phys.num-lattice-sites"].as<unsigned int>() ) {
      throw logic_error(
        "too many electrons (num-electrons > 2 * num-lattice-sites) "
      );
    }

    if ( vm["sim.update-hop-maxdistance"].as<unsigned int>() == 0 ) {
      throw logic_error(
        "electronic configuration updates need at least nearest neighbor hopping"
      );
    }

    if ( vm["sim.update-hop-maxdistance"].as<unsigned int>() > 3 ) {
      throw logic_error(
        "electronic configuration updates with hopping > 3rd nearest neighbors"
        "are not supported"
      );
    }

    if ( vm["phys.lattice"].as<lattice_t>() == LATTICE_2DSQUARE &&
         !is_perfect_square( vm["phys.num-lattice-sites"].as<unsigned int>() ) ) {
      throw logic_error( "the number of lattice sites must be a perfect square" );
    }

    // TODO: minimum lattice size checks (Robert Rueger, 2012-11-17 22:57)

  } catch ( const logic_error& e ) {
    if ( is_master ) {
      cerr << "Logical error in physical parameters: " << e.what() << endl;
    }
    exit( 1 );
  }

  if ( is_master ) {
    cout << endl;
  }

  return vm;
}


istream& operator>>( std::istream& in, lattice_t& lat )
{
  string token;
  in >> token;
  if ( token == "1dchain" ) {
    lat = LATTICE_1DCHAIN;
  } else if ( token == "2dsquare" ) {
    lat = LATTICE_2DSQUARE;
  } else {
    throw po::validation_error( po::validation_error::invalid_option_value );
  }
  return in;
}
