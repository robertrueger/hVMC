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
#include <chrono>
#include <vector>
#include <utility>

#include <eigen3/Eigen/Core>

#include "detwf.hpp"
#include "jastrow.hpp"
#include "fpctrl.hpp"
#include "lattice.hpp"
#include "lattice_1dchain.hpp"
#include "lattice_2dsquare.hpp"

using namespace std;
namespace po = boost::program_options;



BasicSimResults simrun_basic(
  const Options& opts, const VariationalParameters& vpar,
  const boost::mpi::communicator& mpiflock )
{
  if ( mpiflock.rank() == 0 ) {
    cout << ":: Preparing simulation objects ..." << endl;
  }
  HubbardModelVMC* model = nullptr;
  simrun_basic_prepare( opts, vpar, model, mpiflock );

  if ( mpiflock.rank() == 0 ) {
    cout << ":: Equilibrating the system" << endl;
  }
  model->equilibrate( opts["sim.num-mcs-equil"].as<unsigned int>() );

  if ( mpiflock.rank() == 0 ) {
    cout << ":: Performing Monte Carlo cycle" << endl;
  }
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
  // calculate the number of bins this thread has to calculate
  const unsigned int mybincount =
    opts["sim.num-bins"].as<unsigned int>() / mpiflock.size() +
    (
      mpiflock.rank() == 0 ?
      opts["sim.num-bins"].as<unsigned int>() % mpiflock.size() :
      0
    );

  if ( mpiflock.rank() == 0 ) {
    cout << "   -> Creating E_local output array" << endl;
  }
  BinnedData<fptype> E_l( mybincount, opts["sim.num-binmcs"].as<unsigned int>() );

  if ( mpiflock.rank() == 0 ) {
    cout << "   -> Running Monte Carlo cycle" << endl;
  }

  // start the stopwatch
  chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

  for ( unsigned int bin = 0; bin < mybincount; ++bin ) {

    // show progress
    unsigned int progress_percent =
      static_cast<unsigned int>(
        floor( 100 * static_cast<double>( bin ) /
               static_cast<double>( mybincount - 1 ) )
      );

    if ( mpiflock.rank() == 0 ) {
      cout << "\r"
           << "      Bin " << bin + 1 << "/" << mybincount
           << " (" << progress_percent << "%)";
      cout.flush();
    }

    for ( unsigned int mcs = 0;
          mcs < opts["sim.num-binmcs"].as<unsigned int>();
          ++mcs ) {

      if ( !( bin == 0 && mcs == 0 ) ) {
        model->mcs();
      }
      E_l[bin][mcs] = model->E_l();
    }
  }
  if ( mpiflock.rank() == 0 ) {
    cout << endl;
  }

  // stop the stopwatch and calculate the elapsed time
  chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
  const double total_time = chrono::duration<double>( t2 - t1 ).count();
  const double time_per_mcs
    = total_time / static_cast<double>
      ( mybincount * opts["sim.num-binmcs"].as<unsigned int>() );
  const double time_per_phop
    = time_per_mcs / static_cast<double>
      ( opts["phys.num-electrons"].as<unsigned int>() );

  if ( mpiflock.rank() == 0 ) {
    cout << endl
         << "      Finished in " << total_time << " sec" << endl
         << "      => " << 1.0 / time_per_mcs << " MCS/sec" << endl
         << "      => " << 1.0 / ( time_per_phop * 1000.0 ) << " PHOPs/millisec"
         << endl;
  }

  vector<FPDevStat> devstat;
  devstat.push_back( model->get_W_devstat() );
  devstat.push_back( model->get_T_devstat() );

  if ( mpiflock.rank() == 0 ) {
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
  }

  return E_l;
}
