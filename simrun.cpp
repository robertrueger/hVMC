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
  model->equilibrate( opts["sim.num-mcs-equil"].as<unsigned int>() );

  cout << ":: Performing Monte Carlo cycle" << endl;
  const BinnedData<fptype>& E_l = simrun_basic_mccycle( opts, model );

  // we don't need the simulation objects anymore ...
  delete model;
  model = nullptr;

  cout << ":: Performing statistical analysis of the E_local data" << endl;
  const BinnedDataStatistics& E_l_stat = run_bindat_statanalysis( E_l );

  cout << ":: Returning the results" << endl;
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
  unsigned int rngseed
    = opts.count( "sim.rng-seed" ) ?
      opts["sim.rng-seed"].as<unsigned int>() :
      chrono::system_clock::now().time_since_epoch().count();
  cout << "   -> MT19937 RNG ( seed = " << rngseed << " )" << endl;
  mt19937 rng( rngseed  );

  // Lattice object
  Lattice* lat;
  cout << "   -> Lattice object" << endl;
  if ( opts["phys.lattice"].as<lattice_t>() == LATTICE_1DCHAIN ) {
    lat = new Lattice1DChain(
      opts["phys.num-lattice-sites"].as<unsigned int>() );
  } else {
    // TODO: implement 2d square lattice (Robert Rueger, 2012-10-31 15:52)
    exit( 1 );
    lat = new Lattice1DChain(
      opts["phys.num-lattice-sites"].as<unsigned int>() );
  }
  // the lattice object on the heap will be destroyed by hmodvmc's destructor!

  // a vector of the hopping matrix elements
  vector<fptype> t(3);
  t[0] = opts["phys.nn-hopping"].as<fptype>();
  t[1] = opts["phys.2nd-nn-hopping"].as<fptype>();
  t[2] = opts["phys.3rd-nn-hopping"].as<fptype>();
  
  // Jastrow factor
  cout << "   -> Jastrow factor" << endl;
  Jastrow v( lat );
  // TODO: don't set by hand ... (Robert Rueger, 2012-11-09 18:39)
  v.set( 0, 0, -1.f);
  v.set( 0, 1, -.5f);
  v.set( 0, 2, -.2f);

  // determinantal part of the wavefunction
  cout << "   -> Determinantal part of the wavefunction" << endl;
  const Eigen::MatrixXfp& M
    = wf_nntb( t, opts["phys.num-electrons"].as<unsigned int>(), lat );

  // the Hubbard model object itself
  cout << "   -> HubbardModelVMC object" << endl;
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
        opts["sim.num-updates-until-recalc"].as<unsigned int>()
      );
}



BinnedData<fptype> simrun_basic_mccycle(
  const Options& opts, HubbardModelVMC* const model )
{
  cout << "   -> creating E_local output array" << endl;
  BinnedData<fptype> E_l(
    opts["sim.num-bins"].as<unsigned int>(),
    opts["sim.num-binmcs"].as<unsigned int>()
  );

  for ( unsigned int bin = 0;
        bin < opts["sim.num-bins"].as<unsigned int>();
        ++bin ) {
    for ( unsigned int mcs = 0;
          mcs < opts["sim.num-binmcs"].as<unsigned int>();
          ++mcs ) {

      if ( !( bin == 0 && mcs == 0 ) ) {
        model->mcs();
      }
      E_l[bin][mcs] = model->E_l();
    }
  }

  return E_l;
}
