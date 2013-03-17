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
#include <fstream>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/set.hpp>
#include <boost/filesystem/path.hpp>

#include "serialization_eigen.hpp"
#include "mccresults.hpp"
#include "mccrun.hpp"
#include "msgtags.hpp"
#include "fptype.hpp"
#include "obs.hpp"
#include "varparam.hpp"

using namespace std;
namespace mpi = boost::mpi;
namespace fs  = boost::filesystem;


void sched_master( const Options& opts, const mpi::communicator& mpicomm )
{
  schedmsg_t schedmsg;

  // prepare the initial variational parameters
  Eigen::VectorXfp vpar = get_initial_varparam( opts );

  // add the observables you want to measure to the set
  set<observables_t> obs;
  obs.insert( OBSERVABLE_E );
  obs.insert( OBSERVABLE_DELTAK );
  obs.insert( OBSERVABLE_DELTAK_DELTAKPRIME );
  obs.insert( OBSERVABLE_DELTAK_E );

  // create datastructures to store the history of E
  // and the variational parameters during the optimization
  vector<Eigen::VectorXfp> vpar_hist;
  vpar_hist.push_back( vpar );

  // open output files for the energy and the variational parameters
  ofstream vpar_hist_file((
    opts["output-dir"].as<fs::path>() / "vpar_hist.txt"
  ).string());

  ofstream E_hist_file((
    opts["output-dir"].as<fs::path>() / "E_hist.txt"
  ).string());

  // optimization settings
  // TODO: make those as options (Robert Rueger, 2013-03-16 12:34)
  unsigned int sr_bins_init = opts["sim.num-bins"].as<unsigned int>();
  fptype       sr_dt_init
    = opts["sim.sr-dt"].as<fptype>();
  unsigned int sr_max_refinements
    = opts["sim.sr-max-refinements"].as<unsigned int>();
  unsigned int sr_averaging_cycles
    = opts["sim.sr-averaging-cycles"].as<unsigned int>();

  // helper variables
  unsigned int sr_bins = sr_bins_init;
  fptype       sr_dt = sr_dt_init;
  unsigned int sr_cycles = 0;
  unsigned int sr_refinements = 0;
  unsigned int sr_cycles_since_refinement = 0;
  unsigned int sr_num_vpar_converged = 0;
  vector<unsigned int> sr_vpar_converged( vpar.size(), false );
  unsigned int sr_fullconv_cycles = 0;

  bool finished = false;
  while ( !finished ) {

    // start the stopwatch
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

    // tell everyone that we want to start a Monte Carlo cycle
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
                                ( opts["sim.num-binmcs"].as<unsigned int>() *
                                  sr_bins );
    cout << ":: Simulation has finished in " << total_time << " sec" << endl;
    cout << "   Total performance = "
         << 1.0 / time_per_mcs << " effMCS/sec" << endl << endl;


    // output the results
    cout << ":: Simulation results" << endl << endl;
    cout << "      E = " << res.E->mean << endl;
    cout << "sigma_E = " << res.E->sigma << endl << endl;

    if ( opts.count( "verbose" ) ) {
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
    vpar += sr_dt * dvpar;

    // save the new variational parameters for future reference
    vpar_hist.push_back( vpar );

    ++sr_cycles;
    ++sr_cycles_since_refinement;
    if ( sr_num_vpar_converged <= vpar.size() &&
         sr_cycles_since_refinement >= 10 ) {
      for ( unsigned int k = 0; k < vpar.size(); ++k ) {
        // check if the variational parameter is still drifting
        if ( !sr_vpar_converged.at( k ) ) { // ... if it is not converged already
          // (it is considered converged if the sign of its change
          // has been fluctuating during the last iterations)
          int k_signsum = 0;
          for ( auto it = vpar_hist.rbegin(); it != vpar_hist.rbegin() + 9; ++it ) {
            k_signsum += ( *it )( k ) - ( *(it+1) )( k )  < 0.f ? -1 : +1;
          }
          if ( abs( k_signsum ) < 4 ) {
            sr_vpar_converged.at( k ) = true;
            ++sr_num_vpar_converged;
          }
        }
      }
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
    }
    cout << " Iteration: " << sr_cycles
         << " (" << sr_cycles_since_refinement  << ")" << endl;
    cout << "Refinement: " << sr_refinements << "/" << sr_max_refinements << endl;
    cout << " ConvVPars: " << sr_num_vpar_converged << "/" << vpar.size() << endl;
    cout << "    Status: ";

    if ( sr_num_vpar_converged < vpar.size() ) {
      // vpar not converged ... keep iterating!
      cout << "iterating" << endl;
    } else {
      // all variational parameters converged!
      if ( sr_refinements < sr_max_refinements ) {
        // still some refinement to do ... refine!
        sr_dt *= 0.25f;
        sr_bins *= 4;
        sr_vpar_converged = vector<unsigned int>( vpar.size(), false );
        ++sr_refinements;
        sr_cycles_since_refinement = 0;
        sr_num_vpar_converged = 0;
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

  // everything done, tell everyone to quit!
  schedmsg = SCHEDMSG_EXIT;
  mpi::broadcast( mpicomm, schedmsg, 0 );

  // calculate the average of the converged vpars
  const Eigen::VectorXfp vpar_avg =
  accumulate(
    vpar_hist.rbegin() + 1, vpar_hist.rbegin() + sr_averaging_cycles,
    vpar_hist.back()
  ) / static_cast<fptype>( sr_averaging_cycles );

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
