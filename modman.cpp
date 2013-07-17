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

#include "modman.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include "fptype.hpp"

using namespace std;



ModelManager::ModelManager(
  const mt19937& rng_init,
  const shared_ptr<Lattice>& lat_init,
  const DeterminantalWavefunction& detwf_init,
  const Jastrow& v_init,
  unsigned int Ne_init,
  unsigned int update_hop_maxdist_init,
  const vector<double>& t_init,
  double U_init,
  double W_deviation_target,
  unsigned int updates_until_W_recalc,
  double T_deviation_target,
  unsigned int updates_until_T_recalc,
  boost::optional<const Eigen::VectorXi&> spindex_occ_init )
  : rng( rng_init ),
    lat( lat_init ), detwf( detwf_init ), v( v_init ),
    update_hop_maxdist( update_hop_maxdist_init ),
    t( t_init ), U( U_init ),
    pconf( lat, Ne_init, rng, spindex_occ_init ),
    proposed_pconf_accepted( false ),
    W( lat.get(), detwf, pconf, W_deviation_target, updates_until_W_recalc ),
    T( lat.get(), v, pconf, T_deviation_target, updates_until_T_recalc )
{
  bool pconf_has_overlap = W.init_and_check();
  if ( spindex_occ_init ) {
    // there was an inital configuration specified ...
    proposed_pconf_accepted = pconf_has_overlap;
#if VERBOSE >= 1
    if ( proposed_pconf_accepted ) {
      cout << "ModelManager::ModelManager() : "
           << "proposed initial configuration accepted!" << endl;
    } else {
      cout << "ModelManager::ModelManager() : "
           << "proposed initial configuration rejected!" << endl;
    }
#endif
  }

  while ( pconf_has_overlap == false ) {

#if VERBOSE >= 1
    cout << "ModelManager::ModelManager() : configuration does not have "
         << "an overlap with the determinantal wavefunction -> "
         << "generating a new configuration" << endl;
#endif

    pconf.distribute_random();
    pconf_has_overlap = W.init_and_check();
  }

  T.init();

#if VERBOSE >= 1
  cout << "ModelManager::ModelManager() : state has sufficient "
       << "overlap! -> initial state selection completed!" << endl;
#endif
}



bool ModelManager::check_proposed_pconf_accepted() const
{
  return proposed_pconf_accepted;
}



void ModelManager::finalize_equilibration()
{
  verify( W.init_and_check() );
  T.init();
}



void ModelManager::mcs()
{
#if VERBOSE >= 2
  cout << "ModelManager::mcs() : starting new Monte Carlo step!" << endl;
#endif

  // perform a number of metropolis steps equal to the number of electrons
  for ( unsigned int s = 0; s < lat->L; ++s ) {
#if VERBOSE >= 2
    cout << "ModelManager::mcs() : Metropolis step = " << s << endl;
#endif
    metstep();
  }
}



bool ModelManager::metstep()
{
  // let the electron configuration propose a random hop
  const ParticleHop& phop = pconf.propose_random_hop( update_hop_maxdist );


  // check if the hop is possible (hopto site must be empty)
  if ( phop.possible == false ) {

    // hop is not possible, rejected!
#if VERBOSE >= 2
    cout << "ModelManager::metstep() : hop impossible!" << endl;
#endif
    return false;

  } else { // hop possible!

    const double R_j
      = std::exp(
          ( phop.l < lat->L ? 1.0 : -1.0 ) *
          (
            T.get()( lat->get_index_from_spindex( phop.l ) )
            - T.get()( lat->get_index_from_spindex( phop.k_pos ) )
          ) + v.onsite() - v( phop.l, phop.k_pos )
        );

    const double R_s = W.get()( phop.l, phop.k );

    const double accept_prob = R_j * R_j * R_s * R_s;

#if VERBOSE >= 2
    cout << "ModelManager::metstep() : hop possible -> "
         << "R_j = " << R_j
         << ", sdwf_ratio = " << R_s
         << ", accept_prob = " << accept_prob << endl;
#endif

    if ( accept_prob >= 1.0 ||
         uniform_real_distribution<double>( 0.0, 1.0 )( rng ) < accept_prob ) {

#if VERBOSE >= 2
      cout << "ModelManager::metstep() : hop accepted!" << endl;
#endif

      pconf.do_hop( phop );

      W.update( phop );
      T.update( phop );

      return true;

    } else { // hop possible but rejected!

#if VERBOSE >= 2
      cout << "ModelManager::metstep() : hop rejected!" << endl;
#endif

      return false;
    }
  }
}



double ModelManager::E_l() const
{
  // calculate expectation value of the T part of H
  double E_l_kin = 0.0;

  // loop over different elektrons k
  for ( unsigned int k = 0; k < pconf.Np; ++k ) {

    const Lattice::spindex k_pos = pconf.get_particlenum_pos()[ k ];
    assert( pconf.get_spindex_occ()( k_pos ) == PARTICLE_OCCUPATION_FULL );

    // loop over different neighbor orders X
    for ( unsigned int X = 1; X <= t.size(); ++X ) {
      if ( t[X - 1] == 0.0 ) {
        continue;
      }

      double sum_Xnn = 0.0;
      lat->get_Xnn( k_pos, X, &k_pos_Xnn );
      assert( k_pos_Xnn.size() != 0 );

      // loop over different neighbours l of order X
      for ( auto l_it = k_pos_Xnn.begin(); l_it != k_pos_Xnn.end(); ++l_it ) {
        if ( pconf.get_spindex_occ()( *l_it ) == PARTICLE_OCCUPATION_EMPTY ) {

          const double R_j
            = std::exp(
                (
                  ( lat->get_spindex_type( k_pos ) == Lattice::spindex_type::up )
                  ? 1.0 : -1.0
                ) *
                (
                  T.get()( lat->get_index_from_spindex( *l_it ) )
                  - T.get()( lat->get_index_from_spindex( k_pos ) )
                ) + v.onsite() - v( *l_it, k_pos )
              );

          sum_Xnn += R_j * W.get()( *l_it, k );

        }
      }
      // reverse sign of spin down part due to particle-hole-transformation
      E_l_kin +=
        ( lat->get_spindex_type( k_pos ) == Lattice::spindex_type::up ? -1.0 : 1.0 )
        * t[X - 1] * sum_Xnn;

    }
  }

  const double E_l_result =
    ( E_l_kin + U * (
        pconf.get_spindex_occ().head( lat->L ).array() *
        ( 1 - pconf.get_spindex_occ().tail( lat->L ).array() )
      ).sum() ) /
    static_cast<double>( lat->L );

#if VERBOSE >= 2
  cout << "ModelManager::E_l() = " << E_l_result << endl;
#endif

  return E_l_result;
}



Eigen::VectorXd ModelManager::Delta_k( unsigned int optimizers ) const
{
  assert( optimizers > 0 && optimizers < 512 );

  Eigen::VectorXd result = Eigen::VectorXd::Zero( 8 + v.get_num_vpar() );

  // ----- first seven variational parameters are from the determinantal part

  if ( optimizers != 256 ) {
    // only do it if we are optimizing any determinantal parameter

    Eigen::ArrayXfp G = Eigen::ArrayXfp::Zero( 2 * lat->L, 2 * lat->L );
    for ( unsigned int k = 0; k < pconf.Np; ++k ) {
      const Lattice::spindex k_pos = pconf.get_particlenum_pos()[ k ];
      G.row( k_pos ) = W.get().col( k );
    }

    for ( unsigned int vpar = 0; vpar < 8; ++vpar ) {
      if ( ( optimizers >> vpar ) % 2 == 1 ) {
        result( vpar ) = ( detwf.A()[vpar].array() * G ).sum();
      }
    }
  }

  // ----- everything except the first 8 are Jastrow parameters

  if ( optimizers >= 256 ) {
    // only do it if we are optimizing the Jastrow

    for ( Lattice::index i = 0; i < lat->L; ++i ) {
      for ( Lattice::index j = i; j < lat->L; ++j ) {

        const Lattice::irridxrel ij_iir = lat->reduce_idxrel( i, j );
        const double dblcount_correction = ( j == i ) ? 0.5 : 1.0;

        if ( ij_iir != lat->get_maxdist_irridxrel() ) {
          result( 8 + v.get_vparnum( ij_iir ) )
          += dblcount_correction *
             ( pconf.get_spindex_occ()( i ) - pconf.get_spindex_occ()( i + lat->L ) ) *
             ( pconf.get_spindex_occ()( j ) - pconf.get_spindex_occ()( j + lat->L ) );
        }
      }
    }
  }

#if VERBOSE >= 1
  cout << "ModelManager::Delta_k() = " << endl << result.transpose() << endl;
#endif

  return result;
}



double ModelManager::dblocc_dens() const
{
  return
    static_cast<double>(
      (
        pconf.get_spindex_occ().head( lat->L ).array() *
        ( 1 - pconf.get_spindex_occ().tail( lat->L ).array() )
      ).sum()
    ) / static_cast<double>( lat->L );
}



Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> ModelManager::n() const
{
  return (
    pconf.get_spindex_occ().head( lat->L ).array() +
    1 - pconf.get_spindex_occ().tail( lat->L).array()
  ).cast<unsigned int>();
}



Eigen::VectorXi ModelManager::s() const
{
  return (
   pconf.get_spindex_occ().head( lat->L ).array() -
   1 + pconf.get_spindex_occ().tail( lat->L).array()
  );
}


Eigen::VectorXi ModelManager::particleconf() const
{
  return pconf.get_spindex_occ();
}



FPDevStat ModelManager::get_W_devstat() const
{
  return W.get_devstat();
}
FPDevStat ModelManager::get_T_devstat() const
{
  return T.get_devstat();
}
