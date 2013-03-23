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

#include <cmath>
#include <algorithm>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/LU>

#ifdef USE_CBLAS
extern "C" {
# include <cblas.h>
}
#endif

using namespace std;



HubbardModelVMC::HubbardModelVMC(
  const mt19937& rng_init,
  const shared_ptr<Lattice>& lat_init,
  const SingleParticleOrbitals& detwf_init,
  const Jastrow& v_init,
  unsigned int N_init,
  unsigned int update_hop_maxdist_init,
  const vector<fptype>& t_init,
  fptype U_init,
  fptype W_deviation_target,
  unsigned int updates_until_W_recalc,
  fptype T_deviation_target_init,
  unsigned int updates_until_T_recalc_init )
  : rng( rng_init ),
    lat( lat_init ), detwf( detwf_init ), v( v_init ),
    update_hop_maxdist( update_hop_maxdist_init ),
    t( t_init ), U( U_init ),
    econf( ElectronConfiguration( lat, N_init, rng ) ),
    W( WMatrix( lat.get(), detwf, econf,
                W_deviation_target, updates_until_W_recalc ) ),
    T( Eigen::VectorXfp( lat->L ) ),
    updates_until_T_recalc( updates_until_T_recalc_init ),
    updates_since_T_recalc( 0 ),
    T_devstat( FPDevStat( T_deviation_target_init ) )
{

  bool enough_overlap;

  do {

    econf.distribute_random();
    enough_overlap = true; // assume until we are proven wrong

#if VERBOSE >= 2
    cout << "HubbardModelVMC::HubbardModelVMC() : checking newly generated "
         << "state for enough overlap" << endl;
#endif

    // check determinantal part
    enough_overlap &= W.init_and_check();
    if ( !enough_overlap ) {
      continue;
    }

    // check Jastrow part
    T = calc_new_T();
    fptype T_avg = T.squaredNorm() / static_cast<fptype>( T.size() );
    enough_overlap &= T_avg < 100.f;
    if ( !enough_overlap ) {
#if VERBOSE >= 2
      cout << "HubbardModelVMC::HubbardModelVMC() : Jastrow ratios "
           << "are to small, inverse measure is: " << T_avg << endl;
#endif
      continue;
    }

#if VERBOSE >= 2
    cout << "HubbardModelVMC::HubbardModelVMC() : state has sufficient "
         << "overlap! -> initial state selection completed!" << endl;
#endif

  } while ( !enough_overlap );

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

    const fptype R_j = T( lat->get_spinup_site( phop.l ) )
                       / T( lat->get_spinup_site( phop.k_pos ) )
                       * v.exp_onsite() / v.exp( phop.l, phop.k_pos );

    const fptype R_s = W( phop.l, phop.k );

    const fptype accept_prob = R_j * R_j * R_s * R_s;

#if VERBOSE >= 2
    cout << "HubbardModelVMC::metstep() : hop possible -> "
         << "R_j = " << R_j
         << ", sdwf_ratio = " << R_s
         << ", accept_prob = " << accept_prob << endl;
#endif

    if ( accept_prob >= 1.f ||
         uniform_real_distribution<fptype>( 0.f, 1.f )( rng ) < accept_prob ) {

#if VERBOSE >= 2
      cout << "HubbardModelVMC::metstep() : hop accepted!" << endl;
#endif

      econf.do_hop( phop );

      W.update( phop );
      perform_T_update( phop );

      return true;

    } else { // hop possible but rejected!

#if VERBOSE >= 2
      cout << "HubbardModelVMC::metstep() : hop rejected!" << endl;
#endif

      return false;
    }
  }
}



void HubbardModelVMC::perform_T_update( const ElectronHop& hop )
{
  if ( updates_since_T_recalc >= updates_until_T_recalc ) {

#if VERBOSE >= 2
    cout << "HubbardModelVMC::perform_T_update() : recalculating T!" << endl;
#endif

    updates_since_T_recalc = 0;

    const Eigen::MatrixXfp& T_approx = calc_qupdated_T( hop );
    T = calc_new_T();

    fptype dev = calc_deviation( T_approx, T );
    T_devstat.add( dev );

#if VERBOSE >= 2
    cout << "HubbardModelVMC::perform_T_update() : recalculated T "
         << "with deviation = " << dev << endl;

    if ( dev > T_devstat.target ) {
      cout << "HubbardModelVMC::perform_T_update() : deviation goal for matrix "
           << "T not met!" << endl
           << "HubbardModelVMC::perform_T_update() : approximate T =" << endl
           << T_approx.transpose() << endl
           << "HubbardModelVMC::perform_T_update() : exact T =" << endl
           << T.transpose() << endl;
    }
#endif

    assert( dev < T_devstat.target );

  } else {

#if VERBOSE >= 2
    cout << "HubbardModelVMC::perform_T_update() : "
         << "performing a quick update of T!" << endl;
#endif

    ++updates_since_T_recalc;

    T = calc_qupdated_T( hop );

#ifndef NDEBUG
    const Eigen::MatrixXfp& T_chk = calc_new_T();
    fptype dev = calc_deviation( T, T_chk );

# if VERBOSE >= 2
    cout << "HubbardModelVMC::perform_T_update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > T_devstat.target ) {
      cout << "HubbardModelVMC::perform_T_update() : deviation goal for matrix "
           << "T not met!" << endl
           << "HubbardModelVMC::perform_T_update() : quickly updated T =" << endl
           << T.transpose() << endl
           << "HubbardModelVMC::perform_T_update() : exact T =" << endl
           << T_chk.transpose() << endl;
    }
# endif
#endif

    assert( dev < T_devstat.target );
  }
}



Eigen::VectorXfp HubbardModelVMC::calc_new_T() const
{
  Eigen::VectorXfp T_new( lat->L );

  for ( unsigned int i = 0; i < lat->L; ++i ) {
    fptype sum = 0.f;
    for ( unsigned int j = 0; j < lat->L; ++j ) {
      sum += v( i, j ) * static_cast<fptype>(
               ( econf.get_site_occ( j ) + econf.get_site_occ( j + lat->L ) ) );
    }
    T_new( i ) = exp( sum );
  }

  return T_new;
}



Eigen::VectorXfp HubbardModelVMC::calc_qupdated_T( const ElectronHop& hop ) const
{
  Eigen::VectorXfp T_prime( lat->L );

  for ( unsigned int i = 0; i < lat->L; ++i ) {
    T_prime( i ) = T( i ) * v.exp( i, lat->get_spinup_site( hop.l ) )
                   / v.exp( i, lat->get_spinup_site( hop.k_pos ) );
  }

  return T_prime;
}



fptype HubbardModelVMC::E_l()
{
  // calculate expectation value of the T part of H
  fptype E_l_kin = 0.f;

  // loop over different elektrons k
  for ( unsigned int k = 0; k < econf.N(); ++k ) {

    const unsigned int k_pos = econf.get_electron_pos( k );
    assert( econf.get_site_occ( k_pos ) == ELECTRON_OCCUPATION_FULL );

    // loop over different neighbor orders X
    for ( unsigned int X = 1; X <= t.size(); ++X ) {
      if ( t[X - 1] == 0.f ) {
        continue;
      }

      fptype sum_Xnn = 0.f;
      lat->get_Xnn( k_pos, X, &k_pos_Xnn );

      // calculate part of R_j that is constant for this X and k
      assert( k_pos_Xnn.size() != 0 );
      const fptype R_j_constXk =
        v.exp_onsite() / v.exp( k_pos_Xnn[0], k_pos )
        / T( lat->get_spinup_site( k_pos ) );
      // (it is possible to do the idxrel reduction only for one of the
      // neighbours as it is guaranteed to be the same for all of them)

      // loop over different neighbours l of order X
      for ( auto l_it = k_pos_Xnn.begin(); l_it != k_pos_Xnn.end(); ++l_it ) {
        if ( econf.get_site_occ( *l_it ) == ELECTRON_OCCUPATION_EMPTY ) {
          const fptype R_j = T( lat->get_spinup_site( *l_it ) ) * R_j_constXk;
          sum_Xnn += R_j * W( *l_it, k );
        }
      }
      E_l_kin -= t[X - 1] * sum_Xnn;

    }
  }

  const fptype E_l_result =
    ( E_l_kin + U * econf.get_num_dblocc() ) /
    static_cast<fptype>( lat->L );

#if VERBOSE >= 2
  cout << "HubbardModelVMC::E_l() = " << E_l_result << endl;
#endif

  return E_l_result;
}



Eigen::VectorXfp HubbardModelVMC::Delta_k() const
{
  Eigen::VectorXfp sum = Eigen::VectorXfp::Zero( v.get_num_vpar() );

  for ( unsigned int i = 0; i < lat->L; ++i ) {
    for ( unsigned int j = i; j < lat->L; ++j ) {

      unsigned int irr_idxrel = lat->reduce_idxrel( i, j );
      fptype dblcount_correction = ( j == i ) ? .5f : 1.f;

      if ( irr_idxrel != lat->irreducible_idxrel_maxdist() ) {
        unsigned int vparnum = v.get_vparnum( irr_idxrel );
        sum( vparnum ) += dblcount_correction *
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



fptype HubbardModelVMC::dblocc_dens() const
{
  return
    static_cast<fptype>( econf.get_num_dblocc() ) /
    static_cast<fptype>( lat->L );
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
  return T_devstat;
}
