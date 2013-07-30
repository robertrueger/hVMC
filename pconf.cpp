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

#include "pconf.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

using namespace std;



ParticleConfiguration::ParticleConfiguration(
  const shared_ptr<Lattice>& lat_init, unsigned int Ne_init,  mt19937& rng_init,
  boost::optional<const Eigen::VectorXi&> spindex_occ_init )
  : lat( lat_init ),
    Ne( Ne_init ), Npu( Ne / 2 ), Npd( lat->L - Ne / 2 ), Np( Npu + Npd ),
    spindex_occ(
      spindex_occ_init ?
      spindex_occ_init.get() :
      Eigen::VectorXi::Zero( 2 * lat->L )
    ),
    particlenum_pos( Np ),
    rng( rng_init )
{
  assert( Ne % 2 == 0 );

  if ( !spindex_occ_init ) {
    // distribute everything randomly if no initial distribution was given
    distribute_random();
  } else {
    // initial distribution was given

    // make sure the given initial distribution was sane:
    assert( spindex_occ.size() == 2 * lat->L );
    assert( spindex_occ.head( lat->L ).sum() == static_cast<int>( Npu ) );
    assert( spindex_occ.tail( lat->L ).sum() == static_cast<int>( Npd ) );
    assert( spindex_occ.sum() == static_cast<int>( Np ) );

    // we need to construct the particlenum_pos
    // (... the distribute_random method would have done it otherwise ...)
    reconstr_particlenum_pos();
  }
}



void ParticleConfiguration::reconstr_particlenum_pos()
{
  particlenum_pos.clear();
  for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
    if ( spindex_occ[l] == PARTICLE_OCCUPATION_FULL ) {
      particlenum_pos.push_back( l );
    }
  }
  assert( particlenum_pos.size() == Np );

#if VERBOSE >= 2
  cout << "ParticleConfiguration::reconstr_particlenum_pos() : particle positions are"
       << endl;
  for ( auto it = particlenum_pos.begin(); it != particlenum_pos.end(); ++it ) {
    cout << *it << " ";
  }
  cout << endl;
#endif

}



void ParticleConfiguration::distribute_random()
{
  // let the lattice generate a new random configuration
  spindex_occ = lat->get_random_spindex_occ( Npu, Npd, rng );

#if VERBOSE >= 2
  cout << "ParticleConfiguration::distribute_random() : new config is" << endl;
  cout << spindex_occ.head( lat->L ).transpose() << endl
       << spindex_occ.tail( lat->L ).transpose() << endl;
#endif

  assert( spindex_occ.head( lat->L ).sum() == static_cast<int>( Npu ) );
  assert( spindex_occ.tail( lat->L ).sum() == static_cast<int>( Npd ) );
  assert( spindex_occ.sum() == static_cast<int>( Np ) );

  reconstr_particlenum_pos();
}



ParticleHop ParticleConfiguration::propose_random_hop(
  unsigned int update_hop_maxdist ) const
{
  // hop the betath particle
  const unsigned int beta
    = uniform_int_distribution<unsigned int>( 0, Np - 1 )( rng );

  // find the position of the betath particle
  const Lattice::spindex k = particlenum_pos[beta];

#if VERBOSE >= 2
  cout << spindex_occ.head( lat->L ).transpose() << endl
       << spindex_occ.tail( lat->L ).transpose() << endl;
  cout << "ParticleConfiguration::propose_random_hop() : proposing to hop "
       << "particle " << beta << " from " << k << " to ";
#endif

  assert( spindex_occ[ k ] == PARTICLE_OCCUPATION_FULL );

  // get nearest neighbors of site k
  lat->get_Xnn( k, 1, &k_1nb );
  if ( update_hop_maxdist >= 2 ) {
    lat->get_Xnn( k, 2, &k_2nb );
    if ( update_hop_maxdist == 3 ) {
      lat->get_Xnn( k, 3, &k_3nb );
    } else {
      assert( k_3nb.size() == 0 );
    }
  } else {
    assert( k_2nb.size() == 0 );
    assert( k_3nb.size() == 0 );
  }

  const unsigned int nb_number
    = uniform_int_distribution<unsigned int>
      ( 0, k_1nb.size() + k_2nb.size() + k_3nb.size() - 1 )( rng );

  Lattice::spindex l;
  if ( nb_number < k_1nb.size() ) {
    l = k_1nb[ nb_number ];
  } else {
    if ( nb_number < k_1nb.size() + k_2nb.size() ) {
      assert( nb_number >= k_1nb.size() );
      assert( k_2nb.size() != 0 );
      l = k_2nb[ nb_number - k_1nb.size() ];
    } else {
      assert( nb_number >= k_1nb.size() + k_2nb.size() );
      assert( k_3nb.size() != 0 );
      l = k_3nb[ nb_number - k_1nb.size() - k_2nb.size() ];
    }
  }

#if VERBOSE >= 2
  cout << "[?: ";
  for ( auto it = k_1nb.begin(); it != k_1nb.end(); ++it ) {
    cout << *it << " ";
  }
  for ( auto it = k_2nb.begin(); it != k_2nb.end(); ++it ) {
    cout << *it << " ";
  }
  for ( auto it = k_3nb.begin(); it != k_3nb.end(); ++it ) {
    cout << *it << " ";
  }
  cout << "\b] ";
  cout << l << " (" << ( spindex_occ[ l ] == PARTICLE_OCCUPATION_FULL ? "im" : "" )
       << "possible" << ")" << endl;
#endif

  assert( lat->get_spindex_type( k ) == lat->get_spindex_type( l ) );
  assert( spindex_occ[ l ] == PARTICLE_OCCUPATION_FULL ||
          spindex_occ[ l ] == PARTICLE_OCCUPATION_EMPTY  );

  return ParticleHop( beta, l, k, spindex_occ[ l ] == PARTICLE_OCCUPATION_EMPTY );
}



void ParticleConfiguration::do_hop( const ParticleHop& hop )
{
  assert( hop.possible );
  assert( spindex_occ[ hop.k ] == PARTICLE_OCCUPATION_FULL );
  assert( spindex_occ[ hop.l ] == PARTICLE_OCCUPATION_EMPTY );
  assert( lat->get_spindex_type( hop.k ) == lat->get_spindex_type( hop.l ) );

  spindex_occ[ hop.k ] = PARTICLE_OCCUPATION_EMPTY;
  spindex_occ[ hop.l ] = PARTICLE_OCCUPATION_FULL;

  particlenum_pos[ hop.beta ] = hop.l;

#if VERBOSE >= 2
  cout << "ParticleConfiguration::do_hop() : hopping particle #"
       << hop.beta << " from " << hop.k << " to " << hop.l << endl;
  cout << spindex_occ.head( lat->L ).transpose() << endl
       << spindex_occ.tail( lat->L ).transpose() << endl;
  cout << "ParticleConfiguration::do_hop() : particle positions are"
       << endl;
  for ( auto it = particlenum_pos.begin(); it != particlenum_pos.end(); ++it ) {
    cout << *it << " ";
  }
  cout << endl;
#endif

  assert( spindex_occ.head( lat->L ).sum() == static_cast<int>( Npu ) );
  assert( spindex_occ.tail( lat->L ).sum() == static_cast<int>( Npd ) );
  assert( spindex_occ.sum() == static_cast<int>( Np ) );
}


const Eigen::VectorXi& ParticleConfiguration::get_spindex_occ() const
{
  return spindex_occ;
}


const std::vector<Lattice::spindex>& ParticleConfiguration::get_particlenum_pos() const
{
  return particlenum_pos;
}
