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

#include "econf.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

using namespace std;



ElectronConfiguration::ElectronConfiguration(
  const shared_ptr<Lattice>& lat_init,
  unsigned int electron_number_init,
  const shared_ptr<mt19937>& rng_init )
  : lat( lat_init ), electron_number( electron_number_init ),
    site_occ(
      Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>::Zero( 2 * lat_init->L )
    ),
    electron_pos( electron_number_init ),
    rng( rng_init )
{
  distribute_random();
}



void ElectronConfiguration::reconstr_electron_pos()
{
  electron_pos.clear();
  for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
    if ( site_occ[l] == ELECTRON_OCCUPATION_FULL ) {
      electron_pos.push_back( l );
    }
  }
  assert( electron_pos.size() == electron_number );

#if VERBOSE >= 1
  cout << "ElectronConfiguration::reconstr_electron_pos() : electron positions are"
       << endl;
  for ( auto it = electron_pos.begin(); it != electron_pos.end(); ++it ) {
    cout << *it << " ";
  }
  cout << endl;
#endif

}



void ElectronConfiguration::distribute_random()
{
  // TODO: assert -> exception (Robert Rueger, 2012-10-23 13:50)
  assert( electron_number % 2 == 0 );

  // clear all sites
  site_occ = Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>::Zero( 2 * lat->L );

  // randomly distribute L/2 electrons per spin direction
  while ( site_occ.head( lat->L ).sum() < electron_number / 2 ) {
    site_occ[ uniform_int_distribution<unsigned int>( 0, lat->L - 1 )( *rng ) ]
      = ELECTRON_OCCUPATION_FULL;
  }
  while ( site_occ.tail( lat->L ).sum() < electron_number / 2 ) {
    site_occ[ uniform_int_distribution<unsigned int>( 0, lat->L - 1 )( *rng )
              + lat->L ] = ELECTRON_OCCUPATION_FULL;
  }

#if VERBOSE >= 1
  cout << "ElectronConfiguration::distribute_random() : new config is" << endl;
  cout << site_occ.head( lat->L ).transpose() << endl
       << site_occ.tail( lat->L ).transpose() << endl;
#endif

  assert( site_occ.head( lat->L ).sum() == electron_number / 2 );
  assert( site_occ.tail( lat->L ).sum() == electron_number / 2 );
  assert( site_occ.sum() == electron_number );

  reconstr_electron_pos();
}



ElectronHop ElectronConfiguration::propose_random_hop(
  unsigned int update_hop_maxdist )
{
  // hop the kth electron
  const unsigned int k
    = uniform_int_distribution<unsigned int>( 0, electron_number - 1 )( *rng );

  // find the position of the xth electron
  const unsigned int k_pos = electron_pos[k];

#if VERBOSE >= 1
  cout << site_occ.head( lat->L ).transpose() << endl
       << site_occ.tail( lat->L ).transpose() << endl;
  cout << "ElectronConfiguration::propose_random_hop() : proposing to hop "
       << "electron " << k << " from " << k_pos << " to ";
#endif

  assert( site_occ[ k_pos ] == ELECTRON_OCCUPATION_FULL );

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
      ( 0,
        k_1nb.size() + k_2nb.size() + k_3nb.size() - 1
      )( *rng );

  unsigned int l;
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

#if VERBOSE >= 1
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
  cout << l << " (" << ( site_occ[ l ] == ELECTRON_OCCUPATION_FULL ? "im" : "" )
       << "possible" << ")" << endl;
#endif

  assert( ( k_pos < lat->L && l < lat->L ) ||
          ( k_pos >= lat->L && k_pos < 2 * lat->L &&
            l >= lat->L && l < 2 * lat->L ) );
  assert( site_occ[ l ] == ELECTRON_OCCUPATION_FULL ||
          site_occ[ l ] == ELECTRON_OCCUPATION_EMPTY  );

  return ElectronHop( k, l, k_pos, site_occ[ l ] == ELECTRON_OCCUPATION_EMPTY );
}



void ElectronConfiguration::do_hop( const ElectronHop& hop )
{
  assert( hop.possible );
  assert( site_occ[ hop.k_pos ] == ELECTRON_OCCUPATION_FULL );
  assert( site_occ[ hop.l ]     == ELECTRON_OCCUPATION_EMPTY );
  assert( ( hop.k_pos < lat->L && hop.l < lat->L ) ||
          ( hop.k_pos >= lat->L && hop.k_pos < 2 * lat->L &&
            hop.l >= lat->L && hop.l < 2 * lat->L ) );

  site_occ[ hop.k_pos ] = ELECTRON_OCCUPATION_EMPTY;
  site_occ[ hop.l ] = ELECTRON_OCCUPATION_FULL;

  electron_pos[ hop.k ] = hop.l;

#if VERBOSE >= 1
  cout << "ElectronConfiguration::do_hop() : hopping electron #"
       << hop.k << " from " << hop.k_pos << " to " << hop.l << endl;
  cout << site_occ.head( lat->L ).transpose() << endl
       << site_occ.tail( lat->L ).transpose() << endl;
  cout << "ElectronConfiguration::do_hop() : electron positions are"
       << endl;
  for ( auto it = electron_pos.begin(); it != electron_pos.end(); ++it ) {
    cout << *it << " ";
  }
  cout << endl;
#endif

  assert( site_occ.head( lat->L ).sum() == electron_number / 2 );
  assert( site_occ.tail( lat->L ).sum() == electron_number / 2 );
  assert( site_occ.sum() == electron_number );
}



unsigned int ElectronConfiguration::get_electron_pos( unsigned int k ) const
{
  return electron_pos[ k ];
}



unsigned int ElectronConfiguration::get_site_occ( unsigned int l ) const
{
  return site_occ[ l ];
}



unsigned int ElectronConfiguration::N() const
{
  return electron_number;
}



unsigned int ElectronConfiguration::get_num_dblocc() const
{
  return ( site_occ.head( lat->L ).array() * site_occ.tail( lat->L ).array() ).sum();
}
