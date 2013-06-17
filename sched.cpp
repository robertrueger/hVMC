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
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <random>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
//#include <eigen3/Eigen/Eigenvalues>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>

#include "analysis.hpp"
#include "macros.h"
#include "serialization_eigen.hpp"
#include "mccresults.hpp"
#include "mccrun.hpp"
#include "msgtags.hpp"
#include "obs.hpp"
#include "varparam.hpp"
#include "mktest.hpp"

using namespace std;
namespace mpi = boost::mpi;
namespace fs  = boost::filesystem;
namespace ar  = boost::archive;

void sched_master_opt( const Options& opts, const mpi::communicator& mpicomm );
void sched_master_sim( const Options& opts, const mpi::communicator& mpicomm );
void sched_master_ana( const Options& opts );



void sched_master( const Options& opts, const mpi::communicator& mpicomm )
{
  if ( opts["calc.mode"].as<optmode_t>() == OPTION_MODE_OPTIMIZATION ) {

    sched_master_opt( opts, mpicomm );

  } else if ( opts["calc.mode"].as<optmode_t>() == OPTION_MODE_SIMULATION ) {

    sched_master_sim( opts, mpicomm );

  } else if ( opts["calc.mode"].as<optmode_t>() == OPTION_MODE_ANALYSIS ) {

    sched_master_ana( opts );

  }

  // everything done, tell everyone to quit!
  schedmsg_t schedmsg = SCHEDMSG_EXIT;
  mpi::broadcast( mpicomm, schedmsg, 0 );
}



void sched_master_opt( const Options& opts, const mpi::communicator& mpicomm )
{
  // prepare the initial variational parameters
  Eigen::VectorXd vpar = get_initial_varparam( opts );

  // add the observables you want to measure to the set
  set<observables_t> obs;
  obs.insert( OBSERVABLE_E );
  obs.insert( OBSERVABLE_DELTAK );
  obs.insert( OBSERVABLE_DELTAK_DELTAKPRIME );
  obs.insert( OBSERVABLE_DELTAK_E );
  obs.insert( OBSERVABLE_PARTICLE_CONFIGURATIONS );

  // vector to store evolution of the variational parameters during the opt.
  vector<Eigen::VectorXd> vpar_hist;
  vpar_hist.push_back( vpar );

  if ( opts.count( "verbose" ) ) {
    // clear folder for the machine readable variational parameter snapshots
    fs::remove_all(       opts["calc.working-dir"].as<fs::path>() / "vpar_hist" );
    fs::create_directory( opts["calc.working-dir"].as<fs::path>() / "vpar_hist" );
    // initial variational parameters -> first snapshot
    {
      ofstream vpar_init_file( (
        opts["calc.working-dir"].as<fs::path>() / "vpar_hist" / "0.dat"
      ).string() );
      ar::text_oarchive vpar_init_archive( vpar_init_file );
      vpar_init_archive << vpar;
    }
  }

  // open human readable output files for the energy and the variational parameters
  ofstream vpar_hist_file( (
    opts["calc.working-dir"].as<fs::path>() / "opt_vpar_hist.txt"
  ).string(), ios::app );

  ofstream E_hist_file( (
    opts["calc.working-dir"].as<fs::path>() / "opt_E_hist.txt"
  ).string(), ios::app );

  // write the initial variational parameters to the evolution file
  vpar_hist_file << "0 " << vpar.transpose() << endl;
  // ... and to stdout if verbose
  if ( opts.count( "verbose" ) ) {
    cout << "initial vpar = " << endl << vpar.transpose() << endl << endl;
  }

  // optimization settings
  const unsigned int sr_bins_init = opts["calc.num-bins"].as<unsigned int>();
  const double       sr_dt_init   = opts["calc.sr-dt"].as<double>();
  const unsigned int sr_max_refinements
    = opts["calc.sr-max-refinements"].as<unsigned int>();
  const unsigned int sr_averaging_cycles
    = opts["calc.sr-averaging-cycles"].as<unsigned int>();
  const double sr_nodrift_threshold = opts["calc.sr-mkthreshold"].as<double>();
  const double sr_Jboost = opts["calc.sr-dt-Jboost"].as<double>();

  // helper variables
  unsigned int sr_bins = sr_bins_init;
  double       sr_dt = sr_dt_init;
  unsigned int sr_cycles = 0;
  unsigned int sr_refinements = 0;
  unsigned int sr_cycles_since_refinement = 0;
  bool         sr_all_converged = false;
  unsigned int sr_fullconv_cycles = 0;
  vector< MannKendall<double> > sr_vpar_mk( vpar.size() );

  // particle configuration cache
  // (configurations from the last SR iteration used to initialize the next one)
  vector<Eigen::VectorXi> pconf_cache;

  // make a random number generator
  unsigned int rngseed;
  if ( opts.count("calc.rng-seed") ) {
    rngseed = opts["calc.rng-seed"].as<unsigned int>();
  } else {
    rngseed = chrono::system_clock::now().time_since_epoch().count();
  }
  mt19937 rng( rngseed );

  bool finished = false;
  while ( !finished ) {

    // start the stopwatch
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

    // tell everyone that we want to start a Monte Carlo cycle
    schedmsg_t schedmsg;
    schedmsg = SCHEDMSG_START_MCC;
    mpi::broadcast( mpicomm, schedmsg, 0 );

    // send the varparams and the set of observables to the slaves
    mpi::broadcast( mpicomm, vpar, 0 );
    mpi::broadcast( mpicomm, obs,  0 );

    // tell the slaves if we have some good initial confs cached
    bool use_initial_pconf = !pconf_cache.empty();
    mpi::broadcast( mpicomm, use_initial_pconf, 0 );

    MCCResults res;
    if ( use_initial_pconf ) {

      // send a random pconf from the cache to each slave
      for ( int slave = 1; slave < mpicomm.size(); ++slave ) {
        mpicomm.send(
          slave,
          MSGTAG_M_S_INITIAL_PARTICLE_CONFIGURATION,
          pconf_cache[
            uniform_int_distribution<unsigned int>( 0, pconf_cache.size() - 1 )( rng )
          ]
        );
      }

      // run master part of the Monte Carlo cycle without initial pconf
      res =
        mccrun_master(
          opts, vpar, sr_bins, obs, mpicomm,
          pconf_cache[
            uniform_int_distribution<unsigned int>( 0, pconf_cache.size() - 1 )( rng )
          ]
       );

    } else {
      // run master part of the Monte Carlo cycle without initial pconf
      res = mccrun_master( opts, vpar, sr_bins, obs, mpicomm );
    }

    // stop the stopwatch and calculate the elapsed time
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    const double total_time = chrono::duration<double>( t2 - t1 ).count();
    const double time_per_mcs = total_time / static_cast<double>
                                ( opts["calc.num-binmcs"].as<unsigned int>() *
                                  sr_bins );
    cout << ":: Simulation has finished in " << total_time << " sec" << endl;
    cout << "   Total performance = "
         << 1.0 / time_per_mcs << " effMCS/sec" << endl << endl;

    // output the results
    if ( opts.count( "verbose" ) ) {
      cout << ":: Simulation results" << endl;
      cout << res;
    }

    // calculate SR matrix and forces
    const Eigen::MatrixXd S =
      res.Deltak_Deltakprime.get() - res.Deltak.get() * res.Deltak->transpose();
    const Eigen::VectorXd f =
      res.Deltak.get() * res.E->mean - res.Deltak_E.get();

    Eigen::VectorXd dvpar =
      ( S + 0.0001 * Eigen::MatrixXd::Identity( S.rows(), S.cols() ) )
      .fullPivLu().solve( f );

    // Jastrow convergence speed boost
    dvpar.tail( dvpar.size() - 8 ) *= sr_Jboost;

    // update variational parameters
    vpar += sr_dt * dvpar;

    // prevent the absolute value of selected vpars from becoming too small
    for ( unsigned int i = 0; i < 8; ++i ) {
      if ( ( opts["calc.vpar-minabs-select"].as<unsigned int>() >> i ) % 2 == 1 &&
           ( opts["calc.optimizers"].as<unsigned int>() >> i ) % 2 == 1 ) {
        if ( vpar( i ) < +opts["calc.vpar-minabs-value"].as<double>()
             && vpar( i ) > +0.0 ) {
          vpar( i ) = +opts["calc.vpar-minabs-value"].as<double>();
        }
        if ( vpar( i ) > -opts["calc.vpar-minabs-value"].as<double>()
             && vpar( i ) < -0.0 ) {
          vpar( i ) = -opts["calc.vpar-minabs-value"].as<double>();
        }
      }
      assert( std::abs( vpar( i ) )
              >= +opts["calc.vpar-minabs-value"].as<double>() );
    }

    // save the new variational parameters for future reference
    vpar_hist.push_back( vpar );
    // ... and add them to the Mann-Kendall-Test
    for ( unsigned int k = 0; k < vpar.size(); ++k ) {
      sr_vpar_mk[k].push_back( vpar[k] );
    }

    // update particle configuration cache
    pconf_cache = res.pconfs.get();

    ++sr_cycles;
    ++sr_cycles_since_refinement;
    unsigned int sr_num_vpar_converged = 0;
    if ( sr_cycles_since_refinement >= 100 ) {
      for ( unsigned int k = 0; k < vpar.size(); ++k ) {
        // remove old data from the Mann-Kendall test
        // (we want to be testing only the last half since the refinement)
        while ( sr_vpar_mk[k].size() > sr_cycles_since_refinement / 2 ) {
          sr_vpar_mk[k].remove_front();
        }
        // perform a Mann-Kendall test on the evolution of the vpar
        if ( sr_vpar_mk[k].test() <= sr_nodrift_threshold ) {
          ++sr_num_vpar_converged;
        }
      }
    }
    if ( sr_num_vpar_converged == vpar.size() ) {
      sr_all_converged = true;
    }

    // take the machine readable snapshot of the vpars
    if ( opts.count( "verbose" ) ) {
      // determine the file name
      stringstream fname;
      fname << sr_cycles << ".dat";
      ofstream vpar_current_file( (
        opts["calc.working-dir"].as<fs::path>() / "vpar_hist" / fname.str()
      ).string() );
      ar::text_oarchive vpar_current_archive( vpar_current_file );
      vpar_current_archive << vpar;
    }

    // output the vpars and the energy to their files
    vpar_hist_file << sr_cycles << " " << vpar.transpose() << endl;
    E_hist_file << sr_cycles << " " << res.E->mean << " " << res.E->sigma << endl;

    cout << ":: Stochastic reconfiguration" << endl
         << endl;
    if ( opts.count( "verbose" ) ) {
      cout << "S = " << endl << S << endl << endl;
      cout << "f = " << endl << f.transpose() << endl << endl;
      cout << "dvpar = " << endl << dvpar.transpose() << endl << endl;
      cout << "vpar' = " << endl << vpar.transpose() << endl << endl;
      if ( sr_cycles_since_refinement >= 100 ) {
        cout << "mktest = " << endl;
        for ( auto it = sr_vpar_mk.begin(); it != sr_vpar_mk.end(); ++it ) {
          cout << it->test() << " ";
        }
        cout << endl << endl;
      }
    }
    cout << " Iteration: " << sr_cycles
         << " (" << sr_cycles_since_refinement  << ")" << endl;
    cout << "Refinement: " << sr_refinements << "/" << sr_max_refinements << endl;
    cout << " ConvVPars: " << sr_num_vpar_converged << "/" << vpar.size() << endl;
    cout << "    Status: ";

    if ( !sr_all_converged ) {
      // vpar not converged ... keep iterating!
      cout << "iterating" << endl;
    } else {
      // all variational parameters converged!
      if ( sr_refinements < sr_max_refinements ) {
        // still some refinement to do ... refine!
        sr_dt *= 0.5;
        sr_bins *= 2;
        ++sr_refinements;
        sr_cycles_since_refinement = 0;
        sr_all_converged = false;
        cout << "refining & iterating" << endl;
      } else {
        // maximum refinement reached and all vpars fully converged
        ++sr_fullconv_cycles;
        if ( sr_fullconv_cycles < sr_averaging_cycles ) {
          cout << "iterating & measuring" << endl;
        } else {
          cout << "complete" << endl;
          finished = true;
        }
      }
    }

    cout << endl;
  }

  // calculate the average of the converged vpars
  const Eigen::VectorXd vpar_avg =
    accumulate(
      vpar_hist.rbegin() + 1, vpar_hist.rbegin() + sr_averaging_cycles,
      vpar_hist.back()
    ) / static_cast<double>( sr_averaging_cycles );

  // print optimization results
  if ( opts.count( "verbose" ) ) {
    cout << "Evolution of the variational parameters:" << endl;
    for ( unsigned int t = 0; t < vpar_hist.size(); ++t ) {
      cout << t << " " << vpar_hist.at( t ).transpose() << endl;
    }
    cout << endl;
  }
  cout << "Final variational parameters:" << endl;
  cout << vpar_avg.transpose() << endl;

  // write the final variational parameters to a file
  ofstream vpar_final_file( (
    opts["calc.working-dir"].as<fs::path>() / "opt_vpar_final.dat"
  ).string() );
  ar::text_oarchive vpar_final_archive( vpar_final_file );
  vpar_final_archive << vpar_avg;
}



void sched_master_sim( const Options& opts, const mpi::communicator& mpicomm )
{
  schedmsg_t schedmsg;

  // prepare the initial variational parameters
  Eigen::VectorXd vpar = get_initial_varparam( opts );

  // dump the used vpars to disk
  {
    ofstream vpar_file( (
      opts["calc.working-dir"].as<fs::path>() / "sim_vpar.dat"
    ).string() );
    ar::text_oarchive vpar_archive( vpar_file );
    vpar_archive << vpar;
  }

  // add the observables you want to measure to the set
  set<observables_t> obs;
  obs.insert(
    opts["calc.observable"].as< std::vector<observables_t> >().begin(),
    opts["calc.observable"].as< std::vector<observables_t> >().end()
  );

  // start the stopwatch
  chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

  // tell everyone that we want to start a Monte Carlo cycle
  schedmsg = SCHEDMSG_START_MCC;
  mpi::broadcast( mpicomm, schedmsg, 0 );

  // send the varparams and the set of observables to the slaves
  mpi::broadcast( mpicomm, vpar, 0 );
  mpi::broadcast( mpicomm, obs,  0 );

  // tell the slaves if we have some good initial confs cached
  bool use_initial_pconf = false;
  mpi::broadcast( mpicomm, use_initial_pconf, 0 );
  // TODO: read from file or something??? (Robert Rueger, 2013-05-26 12:25)

  // run master part of the Monte Carlo cycle
  const MCCResults& res =
    mccrun_master(
      opts, vpar,
      opts["calc.num-bins"].as<unsigned int>(),
      obs, mpicomm
    );

  // stop the stopwatch and calculate the elapsed time
  chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
  const double total_time = chrono::duration<double>( t2 - t1 ).count();
  const double time_per_mcs = total_time / static_cast<double>
                              ( opts["calc.num-binmcs"].as<unsigned int>() *
                                opts["calc.num-bins"].as<unsigned int>() );
  cout << ":: Simulation has finished in " << total_time << " sec" << endl;
  cout << "   Total performance = "
       << 1.0 / time_per_mcs << " effMCS/sec" << endl << endl;

  // output the results
  if ( opts.count( "verbose" ) ) {
    cout << ":: Simulation results" << endl;
    cout << res;
  }

  // write simulation results to files
  res.write_to_files( opts["calc.working-dir"].as<fs::path>() );
}



void sched_master_ana( const Options& opts )
{
  // figure out the analysis we want to perform ...
  const vector<analysis_t>& anamod_v
    = opts["calc.analysis"].as< vector<analysis_t> >();
  const set<analysis_t> anamod( anamod_v.begin(), anamod_v.end() );

  // ... and hand over control to the analysis modules
  if ( anamod.count( ANALYSIS_STATIC_STRUCTURE_FACTOR ) ) {
    analysis_static_structure_factor(
      opts, opts["calc.working-dir"].as<fs::path>()
    );
  }
}



void sched_slave( const Options& opts, const mpi::communicator& mpicomm )
{
  schedmsg_t schedmsg;
  do {
    mpi::broadcast( mpicomm, schedmsg, 0 );

    if ( schedmsg == SCHEDMSG_START_MCC ) {

      // get variational parameters and set of observables from master
      Eigen::VectorXd vpar;
      mpi::broadcast( mpicomm, vpar, 0 );
      set<observables_t> obs;
      mpi::broadcast( mpicomm, obs,  0 );

      // check if master wants us to initialize the system in a specific pconf
      bool use_initial_pconf;
      mpi::broadcast( mpicomm, use_initial_pconf, 0 );
      if ( use_initial_pconf ) {
        Eigen::VectorXi initial_pconf;
        mpicomm.recv( 0, MSGTAG_M_S_INITIAL_PARTICLE_CONFIGURATION, initial_pconf );

        // run slave part of the Monte Carlo cycle with initial pconf
        mccrun_slave( opts, vpar, obs, mpicomm, initial_pconf );
      } else {
        // run slave part of the Monte Carlo cycle without initial pconf
        mccrun_slave( opts, vpar, obs, mpicomm );
      }
    }

  } while ( schedmsg != SCHEDMSG_EXIT );
}
