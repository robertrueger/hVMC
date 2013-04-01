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

#include "hmodvmc.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

using namespace std;



HubbardModelVMC::HubbardModelVMC(
  const mt19937& rng_init,
  const shared_ptr<Lattice>& lat_init,
  const SingleParticleOrbitals& detwf_init,
  const Jastrow& v_init,
  unsigned int N_init,
  unsigned int update_hop_maxdist_init,
  const vector<double>& t_init,
  double U_init,
  double W_deviation_target,
  unsigned int updates_until_W_recalc,
  double T_deviation_target,
  unsigned int updates_until_T_recalc )
  : rng( rng_init ),
    lat( lat_init ), detwf( detwf_init ), v( v_init ),
    update_hop_maxdist( update_hop_maxdist_init ),
    t( t_init ), U( U_init ),
    econf( ElectronConfiguration( lat, N_init, rng ) ),
    W( WMatrix( lat.get(), detwf, econf,
                W_deviation_target, updates_until_W_recalc ) ),
    T( TVector( lat.get(), v, econf,
                T_deviation_target, updates_until_T_recalc ) )
{
  while ( true ) {

    econf.distribute_random();

#if VERBOSE >= 2
    cout << "HubbardModelVMC::HubbardModelVMC() : checking newly generated "
         << "state for enough overlap" << endl;
#endif

    // check determinantal part for enough overlap
    if ( W.init_and_check() == false ) {
      continue;
    } else {
      break;
    }
  }

  T.init();

#if VERBOSE >= 2
  cout << "HubbardModelVMC::HubbardModelVMC() : state has sufficient "
       << "overlap! -> initial state selection completed!" << endl;
#endif
}



void HubbardModelVMC::mcs()
{
#if VERBOSE >= 2
  cout << "HubbardModelVMC::mcs() : starting new Monte Carlo step!" << endl;
#endif

  // perform a number of metropolis steps equal to the number of electrons
  for ( unsigned int s = 0; s < lat->L; ++s ) {
#if VERBOSE >= 2
    cout << "HubbardModelVMC::mcs() : Metropolis step = " << s << endl;
#endif
    metstep();
  }
}



bool HubbardModelVMC::metstep()
{
  // let the electron configuration propose a random hop
  const ElectronHop& phop = econf.propose_random_hop( update_hop_maxdist );


  // check if the hop is possible (hopto site must be empty)
  if ( phop.possible == false ) {

    // hop is not possible, rejected!
#if VERBOSE >= 2
    cout << "HubbardModelVMC::metstep() : hop impossible!" << endl;
#endif
    return false;

  } else { // hop possible!

    const double R_j = T( lat->get_spinup_site( phop.l ) )
                       / T( lat->get_spinup_site( phop.k_pos ) )
                       * v.exp_onsite() / v.exp( phop.l, phop.k_pos );

    const double R_s = W( phop.l, phop.k );

    const double accept_prob = R_j * R_j * R_s * R_s;

#if VERBOSE >= 2
    cout << "HubbardModelVMC::metstep() : hop possible -> "
         << "R_j = " << R_j
         << ", sdwf_ratio = " << R_s
         << ", accept_prob = " << accept_prob << endl;
#endif

    if ( accept_prob >= 1.0 ||
         uniform_real_distribution<double>( 0.0, 1.0 )( rng ) < accept_prob ) {

#if VERBOSE >= 2
      cout << "HubbardModelVMC::metstep() : hop accepted!" << endl;
#endif

      econf.do_hop( phop );

      W.update( phop );
      T.update( phop );

      return true;

    } else { // hop possible but rejected!

#if VERBOSE >= 2
      cout << "HubbardModelVMC::metstep() : hop rejected!" << endl;
#endif

      return false;
    }
  }
}



double HubbardModelVMC::E_l() const
{
  // calculate expectation value of the T part of H
  double E_l_kin = 0.0;

  // loop over different elektrons k
  for ( unsigned int k = 0; k < econf.N(); ++k ) {

    const unsigned int k_pos = econf.get_electron_pos( k );
    assert( econf.get_site_occ( k_pos ) == ELECTRON_OCCUPATION_FULL );

    // loop over different neighbor orders X
    for ( unsigned int X = 1; X <= t.size(); ++X ) {
      if ( t[X - 1] == 0.0 ) {
        continue;
      }

      double sum_Xnn = 0.0;
      lat->get_Xnn( k_pos, X, &k_pos_Xnn );

      // calculate part of R_j that is constant for this X and k
      assert( k_pos_Xnn.size() != 0 );
      const double R_j_constXk =
        v.exp_onsite() / v.exp( k_pos_Xnn[0], k_pos )
        / T( lat->get_spinup_site( k_pos ) );
      // (it is possible to do the idxrel reduction only for one of the
      // neighbours as it is guaranteed to be the same for all of them)

      // loop over different neighbours l of order X
      for ( auto l_it = k_pos_Xnn.begin(); l_it != k_pos_Xnn.end(); ++l_it ) {
        if ( econf.get_site_occ( *l_it ) == ELECTRON_OCCUPATION_EMPTY ) {
          const double R_j = T( lat->get_spinup_site( *l_it ) ) * R_j_constXk;
          sum_Xnn += R_j * W( *l_it, k );
        }
      }
      E_l_kin -= t[X - 1] * sum_Xnn;

    }
  }

  const double E_l_result =
    ( E_l_kin + U * econf.get_num_dblocc() ) /
    static_cast<double>( lat->L );

#if VERBOSE >= 2
  cout << "HubbardModelVMC::E_l() = " << E_l_result << endl;
#endif

  return E_l_result;
}



Eigen::VectorXd HubbardModelVMC::Delta_k() const
{
  Eigen::VectorXd sum = Eigen::VectorXd::Zero( v.get_num_vpar() );

  for ( unsigned int i = 0; i < lat->L; ++i ) {
    for ( unsigned int j = i; j < lat->L; ++j ) {

      unsigned int irr_idxrel = lat->reduce_idxrel( i, j );
      double dblcount_correction = ( j == i ) ? 0.5 : 1.0;

      if ( irr_idxrel != lat->irreducible_idxrel_maxdist() ) {
        unsigned int vparnum = v.get_vparnum( irr_idxrel );
        sum( vparnum )
        += dblcount_correction *
           ( econf.get_site_occ( i ) + econf.get_site_occ( i + lat->L ) ) *
           ( econf.get_site_occ( j ) + econf.get_site_occ( j + lat->L ) );
      }
    }
  }

#if VERBOSE >= 1
  cout << "HubbardModelVMC::Delta_k() = " << endl << sum.transpose() << endl;
#endif

  return sum;
}



double HubbardModelVMC::dblocc_dens() const
{
  return
    static_cast<double>( econf.get_num_dblocc() ) /
    static_cast<double>( lat->L );
}



Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> HubbardModelVMC::n() const
{
  return econf.n();
}



FPDevStat HubbardModelVMC::get_W_devstat() const
{
  return W.get_devstat();
}
FPDevStat HubbardModelVMC::get_T_devstat() const
{
  return T.get_devstat();
}
