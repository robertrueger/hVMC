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
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

#include <boost/chrono.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/set.hpp>
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
namespace chrono = boost::chrono;
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

  // create datastructures to store the history of E
  // and the variational parameters during the optimization
  vector<Eigen::VectorXd> vpar_hist;
  vpar_hist.push_back( vpar );

  // open output files for the energy and the variational parameters
  ofstream vpar_hist_file( (
    opts["calc.working-dir"].as<fs::path>() / "opt_vpar_hist.txt"
  ).string(), ios::app );

  ofstream E_hist_file( (
    opts["calc.working-dir"].as<fs::path>() / "opt_E_hist.txt"
  ).string(), ios::app );

  // optimization settings
  const unsigned int sr_bins_init = opts["calc.num-bins"].as<unsigned int>();
  const double       sr_dt_init   = opts["calc.sr-dt"].as<double>();
  const unsigned int sr_max_refinements
    = opts["calc.sr-max-refinements"].as<unsigned int>();
  const unsigned int sr_averaging_cycles
    = opts["calc.sr-averaging-cycles"].as<unsigned int>();
  const double sr_nodrift_threshold = opts["calc.sr-mkthreshold"].as<double>();

  // helper variables
  unsigned int sr_bins = sr_bins_init;
  double       sr_dt = sr_dt_init;
  unsigned int sr_cycles = 0;
  unsigned int sr_refinements = 0;
  unsigned int sr_cycles_since_refinement = 0;
  bool         sr_all_converged = false;
  unsigned int sr_fullconv_cycles = 0;
  vector< MannKendall<double> > sr_vpar_mk( vpar.size() );

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

    // run master part of the Monte Carlo cycle
    const MCCResults& res = mccrun_master( opts, vpar, sr_bins, obs, mpicomm );

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
    const Eigen::VectorXd dvpar =
      ( S + 0.01 * Eigen::MatrixXd::Identity( S.rows(), S.cols() ) )
      .fullPivLu().solve( f );

    // update variational parameters
    vpar += sr_dt * dvpar;

    // save the new variational parameters for future reference
    vpar_hist.push_back( vpar );
    // ... and add them to the Mann-Kendall-Test
    for ( unsigned int k = 0; k < vpar.size(); ++k ) {
      sr_vpar_mk[k].push_back( vpar[k] );
    }

    ++sr_cycles;
    ++sr_cycles_since_refinement;
    unsigned int sr_num_vpar_converged = 0;
    if ( sr_cycles_since_refinement >= 20 ) {
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
      if ( sr_cycles_since_refinement >= 20 ) {
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

  // run master part of the Monte Carlo cycle
  const MCCResults& res = mccrun_master(
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

  // write simulation results to human readable files
  res.write_to_files( opts["calc.working-dir"].as<fs::path>() );

  // write results to a machine readable file
  ofstream res_file( (
    opts["calc.working-dir"].as<fs::path>() / "sim_res.dat"
  ).string() );
  ar::text_oarchive res_archive( res_file );
  res_archive << res;
}



void sched_master_ana( const Options& opts )
{
  // read the old MCCResults from disk
  if ( fs::exists(
         opts["calc.working-dir"].as<fs::path>() / "sim_res.dat"
       ) == false ) {
    cout << "ERROR: no simulation result file found" << endl;
    return;
  }
  ifstream res_file( (
    opts["calc.working-dir"].as<fs::path>() / "sim_res.dat"
  ).string() );
  ar::text_iarchive res_archive( res_file );
  MCCResults res;
  res_archive >> res;

  // figure out the analysis we want to perform ...
  const vector<analysis_t>& anamod_v
    = opts["calc.analysis"].as< vector<analysis_t> >();
  const set<analysis_t> anamod( anamod_v.begin(), anamod_v.end() );

  // ... and finally perform the analysis
  if ( anamod.count( ANALYSIS_STATIC_STRUCTURE_FACTOR ) ) {
    analysis_static_structure_factor( opts, res );
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

      // run slave part of the Monte Carlo cycle
      mccrun_slave( opts, vpar, obs, mpicomm );

    }

  } while ( schedmsg != SCHEDMSG_EXIT );
}
