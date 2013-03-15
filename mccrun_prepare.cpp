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

#include "mccrun_prepare.hpp"

#include <vector>

#include "lattice_1dchain.hpp"
#include "lattice_2dsquare.hpp"
#include "obs_all.hpp"

using namespace std;
namespace mpi = boost::mpi;


HubbardModelVMC prepare_model(
  const Options& opts, const Eigen::VectorXfp& vpar,
  const mpi::communicator& mpicomm )
{
  shared_ptr<mt19937> rng = prepare_rng( opts, mpicomm );
  shared_ptr<Lattice> lat = prepare_lattice( opts );

  vector<fptype> t(3);
  t[0] = opts["phys.nn-hopping"].as<fptype>();
  t[1] = opts["phys.2nd-nn-hopping"].as<fptype>();
  t[2] = opts["phys.3rd-nn-hopping"].as<fptype>();

  SingleParticleOrbitals detwf = prepare_detwf( lat, opts );
  if ( detwf.E( detwf.M.cols() ) - detwf.E( detwf.M.cols() - 1 ) < 0.00001 ) {
    if ( mpicomm.rank() == 0 ) {
      cout << endl;
      cout << "   WARNING: Open shell detected!" << endl;
      cout << "     E_fermi = " << detwf.E( detwf.M.cols() - 1 ) << endl;
      cout << "     Orbital below = " << detwf.E( detwf.M.cols() - 2 ) << endl;
      cout << "     Orbital above = " << detwf.E( detwf.M.cols() ) << endl;
      cout << endl;
    }
  }

  Jastrow v( lat, vpar );
  // TODO: REMOVE (Robert Rueger, 2013-01-17 13:31)
  // set some short range terms
  v.set( 0, 0, -1.f );
  v.set( 0, 1, -0.25f );

  return HubbardModelVMC(
    rng,
    lat,
    detwf,
    v,
    opts["phys.num-electrons"].as<unsigned int>(),
    opts["sim.update-hop-maxdistance"].as<unsigned int>(),
    t,
    opts["phys.onsite-energy"].as<fptype>(),
    opts["fpctrl.W-deviation-target"].as<fptype>(),
    opts["fpctrl.W-updates-until-recalc"].as<unsigned int>(),
    opts["fpctrl.T-deviation-target"].as<fptype>(),
    opts["fpctrl.T-updates-until-recalc"].as<unsigned int>()
  );
}


shared_ptr<mt19937> prepare_rng(
  const Options& opts, const mpi::communicator& mpicomm )
{
  unsigned int rngseed = opts["sim.rng-seed"].as<unsigned int>();
  rngseed += rngseed / ( mpicomm.rank() + 1 );
  return make_shared<mt19937>( rngseed );
}


shared_ptr<Lattice> prepare_lattice( const Options& opts )
{
  if ( opts["phys.lattice"].as<lattice_t>() == LATTICE_1DCHAIN ) {
    return make_shared<Lattice1DChain>(
             opts["phys.num-lattice-sites"].as<unsigned int>()
           );
  } else {
    return make_shared<Lattice2DSquare>(
             opts["phys.num-lattice-sites"].as<unsigned int>()
           );
  }
}


SingleParticleOrbitals prepare_detwf(
  const shared_ptr<Lattice>& lat, const Options& opts )
{
  vector<fptype> t(3);
  t[0] = opts["phys.nn-hopping"].as<fptype>();
  t[1] = opts["phys.2nd-nn-hopping"].as<fptype>();
  t[2] = opts["phys.3rd-nn-hopping"].as<fptype>();

  return
    wf_tight_binding( t, opts["phys.num-electrons"].as<unsigned int>(), lat );
}


vector< unique_ptr<Observable> > prepare_obscalcs( const set<observables_t>& obs )
{
  vector< unique_ptr<Observable> > obscalc;

  if ( obs.count( OBSERVABLE_E ) ) {
    obscalc.push_back( unique_ptr<Observable>( new ObservableEnergy() ) );
  }

  if ( obs.count( OBSERVABLE_DELTAK ) ) {
    obscalc.push_back( unique_ptr<Observable>( new ObservableDeltaK() ) );
  }

  if ( obs.count( OBSERVABLE_DELTAK_DELTAKPRIME ) ) {
    obscalc.push_back(
      unique_ptr<Observable>( new ObservableDeltaKDeltaKPrime() )
    );
  }

  if ( obs.count( OBSERVABLE_DELTAK_E ) ) {
    obscalc.push_back(
      unique_ptr<Observable>( new ObservableDeltaKEnergy() )
    );
  }

  return obscalc;
}
