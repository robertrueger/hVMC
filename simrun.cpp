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

#include "simrun.hpp"

#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <functional>

#include <boost/chrono.hpp>
#include <boost/optional.hpp>
#include <boost/mpi/status.hpp>
#include <boost/mpi/request.hpp>

#include <eigen3/Eigen/Core>

#include "detwf.hpp"
#include "jastrow.hpp"
#include "fpctrl.hpp"
#include "lattice.hpp"
#include "lattice_1dchain.hpp"
#include "lattice_2dsquare.hpp"

using namespace std;
namespace po = boost::program_options;
namespace chrono = boost::chrono;
namespace mpi = boost::mpi;



BasicSimResults simrun_basic(
  const Options& opts, const VariationalParameters& vpar,
  const boost::mpi::communicator& mpiflock )
{
  if ( mpiflock.rank() == 0 ) {
    cout << ":: Preparing simulation objects ..." << endl;
  }
  HubbardModelVMC* model = nullptr;
  simrun_basic_prepare( opts, vpar, model, mpiflock );

  const BinnedData<fptype>& E_l = simrun_basic_mccycle( opts, model, mpiflock );

  // we don't need the simulation objects anymore ...
  delete model;
  model = nullptr;

  if ( mpiflock.rank() == 0 ) {
    cout << ":: Performing statistical analysis of the E_local data" << endl;
  }
  const BinnedDataStatistics& E_l_stat = run_bindat_statanalysis( E_l );

  BasicSimResults res;
  res.E       = E_l_stat.mean;
  res.sigma_E = E_l_stat.sigma_mean;
  res.var_E_l = E_l_stat.variance;
  res.success = true;
  return res;
}



void simrun_basic_prepare(
  const Options& opts, const VariationalParameters& vpar,
  HubbardModelVMC*& model, const boost::mpi::communicator& mpiflock )
{
  // Mersenne Twister random number generator
  unsigned int rngseed
    = opts.count( "sim.rng-seed" ) ?
      opts["sim.rng-seed"].as<unsigned int>() :
      chrono::system_clock::now().time_since_epoch().count();
  rngseed += rngseed / ( mpiflock.rank() + 1 );
  if ( mpiflock.rank() == 0 ) {
    cout << "   -> MT19937 RNG ( seed = " << rngseed << " )" << endl;
  }
  mt19937 rng( rngseed  );

  // Lattice object
  Lattice* lat;
  if ( mpiflock.rank() == 0 ) {
    cout << "   -> Lattice object" << endl;
  }
  if ( opts["phys.lattice"].as<lattice_t>() == LATTICE_1DCHAIN ) {
    lat = new Lattice1DChain(
      opts["phys.num-lattice-sites"].as<unsigned int>() );
  } else {
    lat = new Lattice2DSquare(
      opts["phys.num-lattice-sites"].as<unsigned int>() );
  }
  // the lattice object on the heap will be destroyed by hmodvmc's destructor!

  // determinantal part of the wavefunction
  if ( mpiflock.rank() == 0 ) {
    cout << "   -> Determinantal part of the wavefunction" << endl;
  }
  const SingleParticleOrbitals& M
    = wf_tight_binding(
        vpar.determinantal,
        opts["phys.num-electrons"].as<unsigned int>(),
        lat
      );

  // check if the system has an open shell
  if ( M.energies( M.orbitals.cols() ) - M.energies( M.orbitals.cols() - 1 )
       < 0.00001 ) {
    if ( mpiflock.rank() == 0 ) {
      cout << endl;
      cout << "      ERROR: Open shell detected!" << endl;
      cout << "      E_fermi = " << M.energies( M.orbitals.cols() - 1 ) << endl;
      cout << "      Orbital below = " << M.energies( M.orbitals.cols() - 2 ) << endl;
      cout << "      Orbital above = " << M.energies( M.orbitals.cols() ) << endl;
      if ( mpiflock.size() > 1 ) {
        mpiflock.abort( 1 );
      }
    }
    exit( 1 );
  }

  // Jastrow factor
  if ( mpiflock.rank() == 0 ) {
    cout << "   -> Jastrow factor" << endl;
  }
  Jastrow v( lat, vpar.jastrow );
  // TODO: REMOVE (Robert Rueger, 2013-01-17 13:31)
  // set some short range terms
  v.set( 0, 0, -1.f );
  v.set( 0, 1, -0.25f );

  // collect hopping matrix elements into a vector
  vector<fptype> t( 3 );
  t[0] = opts["phys.nn-hopping"].as<fptype>();
  t[1] = opts["phys.2nd-nn-hopping"].as<fptype>();
  t[2] = opts["phys.3rd-nn-hopping"].as<fptype>();

  // the Hubbard model object itself
  if ( mpiflock.rank() == 0 ) {
    cout << "   -> HubbardModelVMC object" << endl;
  }
  model =
    new HubbardModelVMC(
    move( rng ),
    lat,
    move( M ),
    move( v ),
    opts["phys.num-electrons"].as<unsigned int>(),
    opts["sim.update-hop-maxdistance"].as<unsigned int>(),
    move( t ),
    opts["phys.onsite-energy"].as<fptype>(),
    opts["fpctrl.W-deviation-target"].as<fptype>(),
    opts["fpctrl.W-updates-until-recalc"].as<unsigned int>(),
    opts["fpctrl.T-deviation-target"].as<fptype>(),
    opts["fpctrl.T-updates-until-recalc"].as<unsigned int>()
  );
}



BinnedData<fptype> simrun_basic_mccycle(
  const Options& opts, HubbardModelVMC* const model,
  const boost::mpi::communicator& mpiflock )
{
  BinnedData<fptype> E_l( 0, opts["sim.num-binmcs"].as<unsigned int>() );

  if ( mpiflock.rank() == 0 ) {
    // ----- thread is leading the flock -----

    unsigned int finished_flock_members = 0;
    unsigned int scheduled_bins = 0;
    unsigned int completed_bins = 0;
    unsigned int enqueued_bins  = opts["sim.num-bins"].as<unsigned int>();

    // define procedure to query the flock for new work requests
    function<void()> mpiquery_work_requests( [&]() {
      while ( boost::optional<mpi::status> status
              = mpiflock.iprobe( mpi::any_source, 0 ) ) {
        // receive the "I'm ready" message and hand out new bins to the source
        mpiflock.recv( status->source(), 0 );
        if ( enqueued_bins > 0 ) {
          mpiflock.isend( status->source(), 1, 1 );
          scheduled_bins += 1;
          enqueued_bins  -= 1;
        } else {
          mpiflock.isend( status->source(), 1, 0 );
          ++finished_flock_members;
        }
      }
    } );

    // define procedure to query the flock for finished work
    function<void()> mpiquery_finished_work( [&]() {
      while ( boost::optional<mpi::status> status
                = mpiflock.iprobe( mpi::any_source, 2 ) ) {
        mpiflock.recv( status->source(), 2 );
        --scheduled_bins;
        ++completed_bins;
      }
    } );

    unsigned int completed_bins_leader = 0;

    cout << ":: Equilibrating the system" << endl;
    for (
      unsigned int mcs = 0;
      mcs < opts["sim.num-mcs-equil"].as<unsigned int>();
      ++mcs
    ) {
      // 1.: take care of the flock
      mpiquery_finished_work();
      mpiquery_work_requests();

      // 2.: leader calculates a Monte Carlo step in his own cycle
      model->mcs();
    }

    cout << ":: Performing Monte Carlo cycle" << endl;
    cout << endl;
    cout << "      Progress:" << endl;

    while ( enqueued_bins > 0 ) {

      cout << '\r' << "        Bin "
           << completed_bins << "/" << opts["sim.num-bins"].as<unsigned int>();
      cout.flush();

      --enqueued_bins;
      ++scheduled_bins;

      E_l.append_empty_bins( 1 );

      for (
        unsigned int mcs = 0;
        mcs < opts["sim.num-binmcs"].as<unsigned int>();
        ++mcs
      ) {
        // 1.: take care of the flock
        mpiquery_finished_work();
        mpiquery_work_requests();

        // 2.: leader calculates a Monte Carlo step in his own cycle
        model->mcs();
        E_l[completed_bins_leader][mcs] = model->E_l();
      }

      --scheduled_bins;
      ++completed_bins_leader;
      ++completed_bins;
    }
    ++finished_flock_members;

    while ( completed_bins != opts["sim.num-bins"].as<unsigned int>() ||
            static_cast<int>( finished_flock_members ) < mpiflock.size() ) {
      if ( boost::optional<mpi::status> status
           = mpiflock.iprobe( mpi::any_source, 2 ) ) {
        mpiflock.recv( status->source(), 2 );
        --scheduled_bins;
        ++completed_bins;

        cout << '\r' << "        Bin "
             << completed_bins << "/" << opts["sim.num-bins"].as<unsigned int>();
        cout.flush();
      }

      if ( boost::optional<mpi::status> status
           = mpiflock.iprobe( mpi::any_source, 0 ) ) {
        // receive the request for more work
        mpiflock.recv( status->source(), 0 );
        // tell him there is no more work
        mpiflock.isend( status->source(), 1, 0 );
        ++finished_flock_members;
      }
    }
    assert( enqueued_bins == 0 );
    assert( scheduled_bins == 0 );


    cout << '\r' << "        Bin "
         << completed_bins << "/" << opts["sim.num-bins"].as<unsigned int>()
         << endl;
    cout.flush();

    // flock leader checks for floating point problems
    // (the flock will probably not have any if the leader doesn't)
    // TODO: check entire flock for fp probs (Robert Rueger, 2013-01-13 16:05)
    vector<FPDevStat> devstat;
    devstat.push_back( model->get_W_devstat() );
    devstat.push_back( model->get_T_devstat() );

    for ( unsigned int i = 0; i < 2; ++i ) {
      cout << endl;
      cout << "      Floating point precision control for ";
      if ( i == 0 ) {
        cout << "matrix W:";
      } else {
        cout << "vector T:";
      }
      cout << endl;
      cout << "        " << devstat[i].recalcs << " recalculations" << endl;
      cout << "        " << devstat[i].misses << " misses ("
           << devstat[i].mag1_misses << " by a factor of 10)" << endl;
      cout << "        " << devstat[i].hits << " hits ("
           << devstat[i].mag1_hits << " better than a factor of 10)" << endl;
    }
    cout << endl;

  } else {
    // ----- thread is following the flock -----

   // equilibrate the system
   for (
      unsigned int mcs = 0;
      mcs < opts["sim.num-mcs-equil"].as<unsigned int>();
      ++mcs
    ) {
      model->mcs();
    }

    unsigned int completed_bins_thisslave = 0;
    bool leader_out_of_work = false;
    unsigned int scheduled_bins_thisslave;
    mpiflock.isend( 0, 0 );
    mpiflock.recv( 0, 1, scheduled_bins_thisslave );
    leader_out_of_work = ( scheduled_bins_thisslave == 0 );

    while ( scheduled_bins_thisslave > 0 ) {

      unsigned int new_scheduled_bins_thisslave;
      mpi::request leader_answer;
      if ( !leader_out_of_work ) {
        // ask the flock leader for more work
        mpiflock.isend( 0, 0 );
        leader_answer = mpiflock.irecv( 0, 1, new_scheduled_bins_thisslave );
      }

      E_l.append_empty_bins( 1 );

      for (
        unsigned int mcs = 0;
        mcs < opts["sim.num-binmcs"].as<unsigned int>();
        ++mcs
      ) {
        model->mcs();
        E_l[completed_bins_thisslave][mcs] = model->E_l();
      }

      // report completion of the work
      mpiflock.isend( 0, 2 );
      ++completed_bins_thisslave;
      --scheduled_bins_thisslave;

      if ( !leader_out_of_work ) {
        // wait for answer from leader concerning the next bin
        leader_answer.wait();
        if ( new_scheduled_bins_thisslave == 1 ) {
          ++scheduled_bins_thisslave;
        } else {
          leader_out_of_work = true;
        }
      }

    }

  }

  return E_l;
}
