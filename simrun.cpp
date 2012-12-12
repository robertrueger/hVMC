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
using namespace std;
namespace po = boost::program_options;



BasicSimResults simrun_basic( const Options& opts )
{
  cout << ":: Preparing simulation objects ..." << endl;
  HubbardModelVMC* model = nullptr;
  simrun_basic_prepare( opts, model );

  cout << ":: Equilibrating the system" << endl;
  model->equilibrate( opts["sim.num-mcs-equil"].as<cl_uint>() );

  cout << ":: Performing Monte Carlo cycle" << endl;
  const BinnedData<cl_fptype>& E_l = simrun_basic_mccycle( opts, model );

  // we don't need the simulation objects anymore ...
  delete model;
  model = nullptr;

  cout << ":: Performing statistical analysis of the E_local data" << endl;
  const BinnedDataStatistics& E_l_stat = run_bindat_statanalysis( E_l );

  BasicSimResults res;
  res.E       = E_l_stat.mean;
  res.sigma_E = E_l_stat.sigma_mean;
  res.var_E_l = E_l_stat.variance;
  res.success = true;
  return res;
}



void simrun_basic_prepare( const Options& opts, HubbardModelVMC*& model )
{
  // Mersenne Twister random number generator
  cl_uint rngseed
    = opts.count( "sim.rng-seed" ) ?
      opts["sim.rng-seed"].as<cl_uint>() :
      chrono::system_clock::now().time_since_epoch().count();
  cout << "   -> MT19937 RNG ( seed = " << rngseed << " )" << endl;
  mt19937 rng( rngseed  );

  // Lattice object
  Lattice* lat;
  cout << "   -> Lattice object" << endl;
  if ( opts["phys.lattice"].as<lattice_t>() == LATTICE_1DCHAIN ) {
    lat = new Lattice1DChain(
      opts["phys.num-lattice-sites"].as<cl_uint>() );
  } else {
    lat = new Lattice2DSquare(
      opts["phys.num-lattice-sites"].as<cl_uint>() );
  }
  // the lattice object on the heap will be destroyed by hmodvmc's destructor!

  // a vector of the hopping matrix elements
  vector<cl_fptype> t( 3 );
  t[0] = opts["phys.nn-hopping"].as<cl_fptype>();
  t[1] = opts["phys.2nd-nn-hopping"].as<cl_fptype>();
  t[2] = opts["phys.3rd-nn-hopping"].as<cl_fptype>();

  // Jastrow factor
  cout << "   -> Jastrow factor" << endl;
  Jastrow v( lat );

  v.randomize( -.075f, 0.f, &rng );
  v.set( 0, 0, -1.f );
  v.set( 0, 1, -.5f );
  v.set( 0, 7, -.2f );
  v.set( 0, 2, -.1f );

  // determinantal part of the wavefunction
  cout << "   -> Determinantal part of the wavefunction" << endl;
  const SingleParticleOrbitals& M
    = wf_tight_binding( t, opts["phys.num-electrons"].as<cl_uint>(), lat );

  // the Hubbard model object itself
  cout << "   -> HubbardModelVMC object" << endl;
  model =
    new HubbardModelVMC(
    move( rng ),
    lat,
    move( M ),
    move( v ),
    opts["phys.num-electrons"].as<cl_uint>(),
    opts["sim.update-hop-maxdistance"].as<cl_uint>(),
    move( t ),
    opts["phys.onsite-energy"].as<cl_fptype>(),
    opts["fpctrl.W-deviation-target"].as<cl_fptype>(),
    opts["fpctrl.W-updates-until-recalc"].as<cl_uint>(),
    opts["fpctrl.T-deviation-target"].as<cl_fptype>(),
    opts["fpctrl.T-updates-until-recalc"].as<cl_uint>()
  );
}



BinnedData<cl_fptype> simrun_basic_mccycle(
  const Options& opts, HubbardModelVMC* const model )
{
  cout << "   -> Creating E_local output array" << endl;
  BinnedData<cl_fptype> E_l(
    opts["sim.num-bins"].as<cl_uint>(),
    opts["sim.num-binmcs"].as<cl_uint>()
  );

  cout << "   -> Running Monte Carlo cycle" << endl;

  for ( cl_uint bin = 0;
        bin < opts["sim.num-bins"].as<cl_uint>();
        ++bin ) {

    // show progress
    cl_uint progress_percent =
      static_cast<cl_uint>(
        floor( 100 * static_cast<double>( bin ) /
               static_cast<double>( opts["sim.num-bins"].as<cl_uint>() - 1 ) )
      );
    cout << "\r"
         << "      Bin " << bin + 1 << "/" << opts["sim.num-bins"].as<cl_uint>()
         << " (" << progress_percent << "%)";
    cout.flush();

    for ( cl_uint mcs = 0;
          mcs < opts["sim.num-binmcs"].as<cl_uint>();
          ++mcs ) {

      if ( !( bin == 0 && mcs == 0 ) ) {
        model->mcs();
      }
      E_l[bin][mcs] = model->E_l();
    }
  }
  cout << endl;

  vector<FPDevStat> devstat;
  devstat.push_back( model->get_W_devstat() );
  devstat.push_back( model->get_T_devstat() );

  for ( cl_uint i = 0; i < 2; ++i ) {
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

  return E_l;
}
