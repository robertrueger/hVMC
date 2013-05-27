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

#include <boost/chrono.hpp>

#include "lattice_1dchain.hpp"
#include "lattice_2dsquare.hpp"
#include "lattice_2dsqr2layer.hpp"
#include "varparam.hpp"
#include "obs_all.hpp"

using namespace std;
namespace mpi = boost::mpi;
namespace chrono = boost::chrono;


HubbardModelVMC prepare_model(
  const Options& opts, const Eigen::VectorXd& vpar,
  const mpi::communicator& mpicomm,
  boost::optional<const Eigen::VectorXi&> spindex_occ_init )
{
  mt19937 rng = prepare_rng( opts, mpicomm );
  shared_ptr<Lattice> lat = prepare_lattice( opts );

  // slice vpar vector into the different components
  const vector<double> t_vpar
    = { opts["phys.nn-hopping"].as<double>(), vpar[0], vpar[1] };
  const vector<double> Delta_vpar( vpar.data() + 2, vpar.data() + 6 );
  const double mu_vpar = vpar[6];

  const DeterminantalWavefunction detwf
    = build_detwf( lat, opts["phys.num-electrons"].as<unsigned int>(),
                   t_vpar,
                   Delta_vpar, opts["phys.pairing-symmetry"].as<optpairsym_t>(),
                   mu_vpar );

  if ( mpicomm.rank() == 0 && detwf.is_openshell() ) {
    cout << endl;
    cout << "   WARNING: Open shell detected!" << endl;
    cout << endl;
  }

  Jastrow v( lat, vpar.tail( vpar.size() - 7 ) );

  vector<double> t(3);
  t[0] = opts["phys.nn-hopping"].as<double>();
  t[1] = opts["phys.2nd-nn-hopping"].as<double>();
  t[2] = opts["phys.3rd-nn-hopping"].as<double>();

  return HubbardModelVMC(
    rng,
    lat,
    detwf,
    v,
    opts["phys.num-electrons"].as<unsigned int>(),
    opts["calc.update-hop-maxdistance"].as<unsigned int>(),
    t,
    opts["phys.onsite-energy"].as<double>(),
    opts["fpctrl.W-deviation-target"].as<double>(),
    opts["fpctrl.W-updates-until-recalc"].as<unsigned int>(),
    opts["fpctrl.T-deviation-target"].as<double>(),
    opts["fpctrl.T-updates-until-recalc"].as<unsigned int>(),
    spindex_occ_init
  );
}


mt19937 prepare_rng(
  const Options& opts, const mpi::communicator& mpicomm )
{
  unsigned int rngseed;
  if ( opts.count("calc.rng-seed") ) {
    rngseed = opts["calc.rng-seed"].as<unsigned int>();
  } else {
    rngseed = chrono::system_clock::now().time_since_epoch().count();
  }

  rngseed += rngseed / ( mpicomm.rank() + 1 );
  return mt19937( rngseed );
}


shared_ptr<Lattice> prepare_lattice( const Options& opts )
{
  if ( opts["phys.lattice"].as<Lattice::type>() == Lattice::type::chain1d ) {
    return make_shared<Lattice1DChain>(
             opts["phys.num-lattice-sites"].as<unsigned int>()
           );

  } else if ( opts["phys.lattice"].as<Lattice::type>()
                == Lattice::type::square2d )  {
    return make_shared<Lattice2DSquare>(
             opts["phys.num-lattice-sites"].as<unsigned int>()
           );

  } else {
    assert(
      opts["phys.lattice"].as<Lattice::type>() == Lattice::type::square2d2layer
    );

    return make_shared<Lattice2DSquare2Layer>(
             opts["phys.num-lattice-sites"].as<unsigned int>()
           );
  }
}


vector< unique_ptr<Observable> > prepare_obscalcs(
  const set<observables_t>& obs, const Options& opts )
{
  vector< unique_ptr<Observable> > obscalc;

  if ( obs.count( OBSERVABLE_E ) ) {
    obscalc.push_back( unique_ptr<Observable>( new ObservableEnergy() ) );
  }

  if ( obs.count( OBSERVABLE_DELTAK ) ) {
    obscalc.push_back(
      unique_ptr<Observable>(
        new ObservableDeltaK(
          get_num_vpars( opts ),
          opts["calc.optimizers"].as<unsigned int>()
        )
      )
    );
  }

  if ( obs.count( OBSERVABLE_DELTAK_DELTAKPRIME ) ) {
    obscalc.push_back(
      unique_ptr<Observable>(
        new ObservableDeltaKDeltaKPrime(
          get_num_vpars( opts ),
          opts["calc.optimizers"].as<unsigned int>()
        )
      )
    );
  }

  if ( obs.count( OBSERVABLE_DELTAK_E ) ) {
    obscalc.push_back(
      unique_ptr<Observable>(
        new ObservableDeltaKEnergy(
          get_num_vpars( opts ),
          opts["calc.optimizers"].as<unsigned int>()
        )
      )
    );
  }

  if ( obs.count( OBSERVABLE_DOUBLE_OCCUPANCY_DENSITY ) ) {
    obscalc.push_back(
      unique_ptr<Observable>( new ObservableDoubleOccupancy() )
    );
  }

  if ( obs.count( OBSERVABLE_DENSITY_DENSITY_CORRELATION ) ) {
    obscalc.push_back(
      unique_ptr<Observable>(
        new ObservableDensityDensityCorrelation(
          opts["phys.num-lattice-sites"].as<unsigned int>()
        )
      )
    );
  }

  if ( obs.count( OBSERVABLE_SPIN_SPIN_CORRELATION ) ) {
    obscalc.push_back(
      unique_ptr<Observable>(
        new ObservableSpinSpinCorrelation(
          opts["phys.num-lattice-sites"].as<unsigned int>()
        )
      )
    );
  }

  if ( obs.count( OBSERVABLE_PARTICLE_CONFIGURATIONS ) ) {
    obscalc.push_back(
      unique_ptr<Observable>( new ObservableParticleConfigurations() )
    );
  }

  return obscalc;
}
