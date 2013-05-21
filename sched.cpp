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
#include <sstream>
#include <algorithm>
#include <iterator>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
//#include <eigen3/Eigen/Eigenvalues>

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

  // vector to store evolution of the variational parameters during the opt.
  vector<Eigen::VectorXd> vpar_hist;
  vpar_hist.push_back( vpar );

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

/*
    // ---- SORELLA METHOD

    // calculate the required change of the variational parameters
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> S_esolver( S + 0.001 * Eigen::MatrixXd::Identity( S.rows(), S.cols() ) );
    assert( S_esolver.info() == Eigen::Success );
    Eigen::VectorXd gamma
      =   ( S_esolver.eigenvectors().adjoint() * f ).array()
        / S_esolver.eigenvalues().array();
    for ( unsigned int i = 0; i < vpar.size(); ++i ) {
      if ( S_esolver.eigenvalues()( i ) / S_esolver.eigenvalues()( vpar.size() - 1 ) < 1.0E-10 ) {
        gamma( i ) = 0;
      }
    }
    Eigen::VectorXd dvpar = S_esolver.eigenvectors() * gamma;
*/


    // ---- TOCCHIO METHOD

    Eigen::VectorXd dvpar =
      ( S + 0.0001 * Eigen::MatrixXd::Identity( S.rows(), S.cols() ) )
      .fullPivLu().solve( f );

/*
    // ---- ATTACCALITE METHOD

    // calculate eigenvalues of S
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> S_esolver( S );
    assert( S_esolver.info() == Eigen::Success );

#if VERBOSE >= 1
    cout << "sched_master_opt() : eigenvalues of S =" << endl
         << S_esolver.eigenvalues().transpose() << endl;
#endif

    // count the number of eigenvalues where the sqrt is smaller than epsilon
    const unsigned int p = ( S_esolver.eigenvalues().array().sqrt() <= 0.001 ).count();

#if VERBOSE >= 1
    cout << "sched_master_opt() : number of eigenvalues below threshold = " << p << endl;
#endif

    // define a reduced matrix S
    Eigen::MatrixXd S_red = S;
    // ... and a vector that keeps track of which variational parameter
    // ends up in which row/column of the S_red matrix and our final solution
    std::vector<unsigned int> red_vpartracker( vpar.size() );
    for ( unsigned int i = 0; i < vpar.size(); ++i ) {
      red_vpartracker.at( i ) = i;
    }

    // loop that removes rows and columns from S; reduces S to S_red
    for ( unsigned int k = 0; k < p; ++k ) {

      // invert S_red
      Eigen::MatrixXd S_red_inv = S_red.fullPivLu().inverse();

#if VERBOSE >= 1
      cout << "sched_master_opt() : S_red_inv =" << endl << S_red_inv << endl;
#endif

      // find index i of the largest diagonal element of S_red_inv
      unsigned int i_max = 0;
      for ( unsigned int i = 1; i < S_red_inv.diagonal().size(); ++i ) {
        if ( S_red_inv.diagonal()( i ) > S_red_inv.diagonal()( i_max ) ) {
          i_max = i;
        }
      }

#if VERBOSE >= 1
      cout << "sched_master_opt() : largest diagonal element of S_red_inv = "
           << S_red_inv.diagonal()( i_max ) << " (at position " << i_max << ")" << endl;
#endif

      // remove the i-th row and column from S_red
      for ( unsigned int j = i_max + 1; j < S_red.rows(); ++j ) {
        S_red.row( j - 1 ) = S_red.row( j );
      }
      for ( unsigned int j = i_max + 1; j < S_red.cols(); ++j ) {
        S_red.col( j - 1 ) = S_red.col( j );
      }
      S_red.conservativeResize( S_red.rows() - 1, S_red.cols() - 1 );

#if VERBOSE >= 1
      cout << "sched_master_opt() : new S_red =" << endl << S_red << endl;
#endif

      // shift the variational parameter tracker
      red_vpartracker.erase( red_vpartracker.begin() + i_max );
      assert( static_cast<unsigned int>( red_vpartracker.size() ) == S_red.rows() );
      assert( static_cast<unsigned int>( red_vpartracker.size() ) == S_red.cols() );

#if VERBOSE >= 1
      cout << "sched_master_opt() : new vpartracker = " << endl;
      for ( auto it = red_vpartracker.begin(); it != red_vpartracker.end(); ++it ) {
        cout << *it << " ";
      }
      cout << endl;
#endif
    }

    // build the reduced force vector
    Eigen::VectorXd f_red( S_red.cols() );
    for ( unsigned int i = 0; i < f_red.size(); ++i ) {
     f_red( i ) = f( red_vpartracker.at( i ) );
    }

#if VERBOSE >= 1
    cout << "sched_master_opt() : f_red =" << endl << f_red.transpose() << endl;
#endif

    // calculate the change of the reduced set of vpars by solving S_red * dv_red = f_red
    const Eigen::VectorXd dvpar_red =
      ( S_red + 0.0001 * Eigen::MatrixXd::Identity( S_red.rows(), S_red.cols() ) )
      .fullPivLu().solve( f_red );

#if VERBOSE >= 1
    cout << "sched_master_opt() : dvpar_red =" << endl << dvpar_red.transpose() << endl;
#endif

    // build change vector from reduced change vector (filling in the gaps with 0)
    Eigen::VectorXd dvpar = Eigen::VectorXd::Zero( vpar.size() );
    for ( unsigned int i = 0; i < dvpar_red.size(); ++i ) {
      dvpar( red_vpartracker.at( i ) ) = dvpar_red( i );
    }

#if VERBOSE >= 1
    cout << "sched_master_opt() : dvpar =" << endl << dvpar.transpose() << endl;
#endif
*/

    // Jastrow convergence speed boost
    dvpar.tail( dvpar.size() - 7 ) *= sr_Jboost;

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

    // take the machine readable snapshot of the vpars
    {
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

      // run slave part of the Monte Carlo cycle
      mccrun_slave( opts, vpar, obs, mpicomm );

    }

  } while ( schedmsg != SCHEDMSG_EXIT );
}
