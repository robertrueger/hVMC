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

#include "mccrun.hpp"

#include <random>
#include <memory>
#include <vector>
#include <functional>
#include <numeric>

#include <boost/mpi/collectives.hpp>

#include "mccrun_prepare.hpp"
#include "hmodvmc.hpp"
#include "obs_all.hpp"
#include "msgtags.hpp"

using namespace std;
namespace mpi = boost::mpi;



MCCResults mccrun_master(
  const Options& opts, const Eigen::VectorXd& vpar, unsigned int num_bins,
  const set<observables_t>& obs, const mpi::communicator& mpicomm )
{
  cout << "========== NEW MONTE CARLO CYCLE ==========" << endl;
  cout << ":: Preparing the simulation" << endl;

  HubbardModelVMC model = prepare_model( opts, vpar, mpicomm );

  vector< unique_ptr<Observable> > obscalc = prepare_obscalcs( obs, opts );
  ObservableCache obscache;

  unsigned int finished_workers = 0;
  unsigned int scheduled_bins = 0;
  unsigned int completed_bins = 0;
  unsigned int enqueued_bins  = num_bins;

  // define procedure to query the slaves for new work requests
  function<void()> mpiquery_work_requests( [&]() {
    while ( boost::optional<mpi::status> status
            = mpicomm.iprobe( mpi::any_source, MSGTAG_S_M_REQUEST_BINS ) ) {
      // receive the request and hand out new bins to the source
      mpicomm.recv( status->source(), MSGTAG_S_M_REQUEST_BINS );
      if ( enqueued_bins > 0 ) {
        mpicomm.send( status->source(), MSGTAG_M_S_DISPATCHED_BINS, 1 );
        scheduled_bins += 1;
        enqueued_bins  -= 1;
      } else {
        mpicomm.send( status->source(), MSGTAG_M_S_DISPATCHED_BINS, 0 );
        ++finished_workers;
      }
    }
  } );

  // define procedure to query the slaves for finished work
  function<void()> mpiquery_finished_work( [&]() {
    while ( boost::optional<mpi::status> status
            = mpicomm.iprobe( mpi::any_source, 2 ) ) {
      mpicomm.recv( status->source(), 2 );
      --scheduled_bins;
      ++completed_bins;
    }
  } );

  cout << ":: Equilibrating the system" << endl;

  for (
    unsigned int mcs = 0;
    mcs < opts["calc.num-mcs-equil"].as<unsigned int>();
    ++mcs ) {
    // take care of the slaves
    mpiquery_finished_work();
    mpiquery_work_requests();

    // perform a Monte Carlo step
    model.mcs();
  }

  unsigned int completed_bins_master = 0;

  cout << ":: Performing Monte Carlo cycle" << endl;
  cout << endl;
  cout << "   Progress:" << endl;

  while ( enqueued_bins > 0 ) {

    cout << '\r' << "     Bin "
         << completed_bins << "/" << num_bins;
    cout.flush();

    --enqueued_bins;
    ++scheduled_bins;

    for (
      unsigned int mcs = 0;
      mcs < opts["calc.num-binmcs"].as<unsigned int>();
      ++mcs ) {
      // take care of the slaves
      mpiquery_finished_work();
      mpiquery_work_requests();

      // perform a Monte Carlo step
      model.mcs();

      // measure observables
      for ( const unique_ptr<Observable>& o : obscalc ) {
        o->measure( model, obscache );
      }
      obscache.clear();
    }

    // tell the observables that a bin has been completed
    for ( const unique_ptr<Observable>& o : obscalc ) {
      o->completebin();
    }

    --scheduled_bins;
    ++completed_bins_master;
    ++completed_bins;
  }
  ++finished_workers;

  while ( completed_bins != num_bins ||
          static_cast<int>( finished_workers ) < mpicomm.size() ) {
    if ( boost::optional<mpi::status> status
         = mpicomm.iprobe( mpi::any_source, MSGTAG_S_M_FINISHED_BINS ) ) {
      mpicomm.recv( status->source(), MSGTAG_S_M_FINISHED_BINS );
      --scheduled_bins;
      ++completed_bins;

      cout << '\r' << "     Bin " << completed_bins << "/" << num_bins;
      cout.flush();
    }

    if ( boost::optional<mpi::status> status
         = mpicomm.iprobe( mpi::any_source, MSGTAG_S_M_REQUEST_BINS ) ) {
      // receive the request for more work
      mpicomm.recv( status->source(), MSGTAG_S_M_REQUEST_BINS );
      // tell him there is no more work
      mpicomm.send( status->source(), MSGTAG_M_S_DISPATCHED_BINS, 0 );
      ++finished_workers;
    }
  }
  assert( enqueued_bins == 0 );
  assert( scheduled_bins == 0 );

  cout << '\r' << "     Bin " << completed_bins << "/" << num_bins << endl;
  cout.flush();

  // check for floating point precision problems

  cout << endl;
  cout << "   Floating point precision control" << endl;

  vector<FPDevStat> W_devstats;
  assert( mpicomm.rank() == 0 );
  mpi::gather( mpicomm, model.get_W_devstat(), W_devstats, 0 );
  FPDevStat W_devstat_combined =
    accumulate(
      W_devstats.begin(), W_devstats.end(),
      FPDevStat( opts["fpctrl.W-deviation-target"].as<double>() )
    );
  cout << "     W: " << W_devstat_combined.recalcs
       << "/" << W_devstat_combined.misses
       << "/" << W_devstat_combined.mag1_misses << endl;

  vector<FPDevStat> T_devstats;
  assert( mpicomm.rank() == 0 );
  mpi::gather( mpicomm, model.get_T_devstat(), T_devstats, 0 );
  FPDevStat T_devstat_combined =
    accumulate(
      T_devstats.begin(), T_devstats.end(),
      FPDevStat( opts["fpctrl.T-deviation-target"].as<double>() )
    );
  cout << "     T: " << T_devstat_combined.recalcs
       << "/" << T_devstat_combined.misses
       << "/" << T_devstat_combined.mag1_misses << endl;

  if ( W_devstat_combined.mag1_misses > 0 ||
       T_devstat_combined.mag1_misses > 0 ) {
    cout << "   Precision targets missed by more than an order of magnitude!" << endl
         << "   WARNING: Your results might be unreliable!!!" << endl << endl;
  } else if ( W_devstat_combined.misses > 0 ||
              T_devstat_combined.misses > 0 ) {
    cout << "   Some precision targets were missed, but your results should be fine."
         << endl << endl;
  } else {
    cout << "   No missed precision targets." << endl << endl;
  }

  // collect results from the slaves and return results the scheduler
  MCCResults results;
  for ( const unique_ptr<Observable>& o : obscalc ) {
    o->collect_and_write_results( mpicomm, results );
  }
  results.success = true;
  return results;
}


void mccrun_slave(
  const Options& opts, const Eigen::VectorXd& vpar,
  const set<observables_t>& obs, const mpi::communicator& mpicomm )
{
  // prepare the simulation

  HubbardModelVMC model = prepare_model( opts, vpar, mpicomm );
  vector< unique_ptr<Observable> > obscalc = prepare_obscalcs( obs, opts );
  ObservableCache obscache;

  // equilibrate the system

  for (
    unsigned int mcs = 0;
    mcs < opts["calc.num-mcs-equil"].as<unsigned int>();
    ++mcs )
  {
    model.mcs();
  }

  // run this slaves part of the Monte Carlo cycle

  unsigned int completed_bins_thisslave = 0;
  bool master_out_of_work = false;
  unsigned int scheduled_bins_thisslave;
  mpicomm.send( 0, MSGTAG_S_M_REQUEST_BINS );
  mpicomm.recv( 0, MSGTAG_M_S_DISPATCHED_BINS, scheduled_bins_thisslave );
  master_out_of_work = ( scheduled_bins_thisslave == 0 );

  while ( scheduled_bins_thisslave > 0 ) {

    unsigned int new_scheduled_bins_thisslave;
    mpi::request master_answer;
    if ( !master_out_of_work ) {
      // ask the master for more work
      mpicomm.send( 0, MSGTAG_S_M_REQUEST_BINS );
      master_answer = mpicomm.irecv(
        0, MSGTAG_M_S_DISPATCHED_BINS,
        new_scheduled_bins_thisslave
      );
    }

    for (
      unsigned int mcs = 0;
      mcs < opts["calc.num-binmcs"].as<unsigned int>();
      ++mcs )
    {
      // perform a Monte Carlo step
      model.mcs();

     // measure observables
      for ( const unique_ptr<Observable>& o : obscalc ) {
        o->measure( model, obscache );
      }
      obscache.clear();
    }

    // tell the observables that a bin has been completed
    for ( const unique_ptr<Observable>& o : obscalc ) {
      o->completebin();
    }

    // report completion of the work
    mpicomm.send( 0, 2 );
    ++completed_bins_thisslave;
    --scheduled_bins_thisslave;

    if ( !master_out_of_work ) {
      // wait for answer from master concerning the next bin
      master_answer.wait();
      if ( new_scheduled_bins_thisslave == 1 ) {
        ++scheduled_bins_thisslave;
      } else {
        master_out_of_work = true;
      }
    }
  }

  // send floating point precision control data to master
  mpi::gather( mpicomm, model.get_W_devstat(), 0 );
  mpi::gather( mpicomm, model.get_T_devstat(), 0 );

  // send observables to master
  for ( const unique_ptr<Observable>& o : obscalc ) {
    o->send_results_to_master( mpicomm );
  }
}
