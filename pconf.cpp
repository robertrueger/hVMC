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
  const shared_ptr<Lattice>& lat_init, unsigned int Ne_init,  mt19937& rng_init )
  : lat( lat_init ),
    Ne( Ne_init ), Npu( Ne / 2 ), Npd( lat->L - Ne / 2 ), Np( Npu + Npd ),
    site_occ( Eigen::VectorXi::Zero( 2 * lat->L ) ),
    particle_pos( Np ),
    rng( rng_init )
{
  assert( Ne % 2 == 0 );
  distribute_random();
}



void ParticleConfiguration::reconstr_particle_pos()
{
  particle_pos.clear();
  for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
    if ( site_occ[l] == PARTICLE_OCCUPATION_FULL ) {
      particle_pos.push_back( l );
    }
  }
  assert( particle_pos.size() == Np );

#if VERBOSE >= 2
  cout << "ParticleConfiguration::reconstr_particle_pos() : particle positions are"
       << endl;
  for ( auto it = particle_pos.begin(); it != particle_pos.end(); ++it ) {
    cout << *it << " ";
  }
  cout << endl;
#endif

}



void ParticleConfiguration::distribute_random()
{
  // let the lattice generate a new random configuration
  site_occ = lat->get_random_site_occ( Npu, Npd, rng );

#if VERBOSE >= 2
  cout << "ParticleConfiguration::distribute_random() : new config is" << endl;
  cout << site_occ.head( lat->L ).transpose() << endl
       << site_occ.tail( lat->L ).transpose() << endl;
#endif

  assert( site_occ.head( lat->L ).sum() == static_cast<int>( Npu ) );
  assert( site_occ.tail( lat->L ).sum() == static_cast<int>( Npd ) );
  assert( site_occ.sum() == static_cast<int>( Np ) );

  reconstr_particle_pos();
}



ParticleHop ParticleConfiguration::propose_random_hop(
  unsigned int update_hop_maxdist ) const
{
  // hop the kth particle
  const unsigned int k
    = uniform_int_distribution<unsigned int>( 0, Np - 1 )( rng );

  // find the position of the kth particle
  const Lattice::spindex k_pos = particle_pos[k];

#if VERBOSE >= 2
  cout << site_occ.head( lat->L ).transpose() << endl
       << site_occ.tail( lat->L ).transpose() << endl;
  cout << "ParticleConfiguration::propose_random_hop() : proposing to hop "
       << "particle " << k << " from " << k_pos << " to ";
#endif

  assert( site_occ[ k_pos ] == PARTICLE_OCCUPATION_FULL );

  // get nearest neighbors of site k_pos
  lat->get_Xnn( k_pos, 1, &k_1nb );
  if ( update_hop_maxdist >= 2 ) {
    lat->get_Xnn( k_pos, 2, &k_2nb );
    if ( update_hop_maxdist == 3 ) {
      lat->get_Xnn( k_pos, 3, &k_3nb );
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
  cout << l << " (" << ( site_occ[ l ] == PARTICLE_OCCUPATION_FULL ? "im" : "" )
       << "possible" << ")" << endl;
#endif

  assert( lat->get_spindex_type( k_pos ) == lat->get_spindex_type( l ) );
  assert( site_occ[ l ] == PARTICLE_OCCUPATION_FULL ||
          site_occ[ l ] == PARTICLE_OCCUPATION_EMPTY  );

  return ParticleHop( k, l, k_pos, site_occ[ l ] == PARTICLE_OCCUPATION_EMPTY );
}



void ParticleConfiguration::do_hop( const ParticleHop& hop )
{
  assert( hop.possible );
  assert( site_occ[ hop.k_pos ] == PARTICLE_OCCUPATION_FULL );
  assert( site_occ[ hop.l ]     == PARTICLE_OCCUPATION_EMPTY );
  assert( lat->get_spindex_type( hop.k_pos ) == lat->get_spindex_type( hop.l ) );

  site_occ[ hop.k_pos ] = PARTICLE_OCCUPATION_EMPTY;
  site_occ[ hop.l ] = PARTICLE_OCCUPATION_FULL;

  particle_pos[ hop.k ] = hop.l;

#if VERBOSE >= 2
  cout << "ParticleConfiguration::do_hop() : hopping particle #"
       << hop.k << " from " << hop.k_pos << " to " << hop.l << endl;
  cout << site_occ.head( lat->L ).transpose() << endl
       << site_occ.tail( lat->L ).transpose() << endl;
  cout << "ParticleConfiguration::do_hop() : particle positions are"
       << endl;
  for ( auto it = particle_pos.begin(); it != particle_pos.end(); ++it ) {
    cout << *it << " ";
  }
  cout << endl;
#endif

  assert( site_occ.head( lat->L ).sum() == static_cast<int>( Npu ) );
  assert( site_occ.tail( lat->L ).sum() == static_cast<int>( Npd ) );
  assert( site_occ.sum() == static_cast<int>( Np ) );
}



Lattice::spindex ParticleConfiguration::get_particle_pos( unsigned int k ) const
{
  return particle_pos[ k ];
}



ParticleOccupation_t ParticleConfiguration::get_site_occ( Lattice::spindex l ) const
{
  return static_cast<ParticleOccupation_t>( site_occ[ l ] );
}



Eigen::VectorXi ParticleConfiguration::npu() const
{
  return site_occ.head( lat->L );
}

Eigen::VectorXi ParticleConfiguration::npd() const
{
  return site_occ.tail( lat->L );
}
