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

#include "econf.hpp"
using namespace std;



ElectronConfiguration::ElectronConfiguration(
  Lattice* const lat_init,
  unsigned int electron_number_init,
  mt19937* rng_init )
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


/*
void ElectronConfiguration::distribute_mixed()
{
  // TODO: assert -> exception (Robert Rueger, 2012-10-23 13:50)
  assert( electron_number % 2 == 0 );

  // clear all sites
  site_occ = Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>::Zero( 2 * lat->L );

  // randomly distribute L/2 electrons per spin direction
  for ( unsigned int l = 0; l < lat->L - 1; l += 2 ) {
    site_occ[l] = ELECTRON_OCCUPATION_FULL;
    site_occ[l + 1 + lat->L] = ELECTRON_OCCUPATION_FULL;
  }

  // TODO: generalize to != 1/2 filling (Robert Rueger, 2012-10-25 17:08)
  assert( electron_number == lat->L );

#if VERBOSE >= 2
  cout << "ElectronConfiguration::distribute_random() : new config is" << endl;
  cout << site_occ.head( lat->L ).transpose() << endl
       << site_occ.tail( lat->L ).transpose() << endl;
#endif

  assert( site_occ.head( lat->L ).sum() == electron_number / 2 );
  assert( site_occ.tail( lat->L ).sum() == electron_number / 2 );
  assert( site_occ.sum() == electron_number );

  reconstr_electron_pos();
}
*/


ElectronHop ElectronConfiguration::propose_random_hop(
  unsigned int update_hop_maxdist ) const
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

  vector<unsigned int> k_nb = lat->get_Xnn( k_pos, 1 );
  // append further away neighbors
  for ( unsigned int X = 2; X <= update_hop_maxdist; ++X ) {
    const vector<unsigned int>& k_Xnn = lat->get_Xnn( k_pos, X );
    k_nb.insert( k_nb.end(), k_Xnn.begin(), k_Xnn.end() );
  }

  const unsigned int l
    = k_nb[ uniform_int_distribution<unsigned int>( 0, k_nb.size() - 1 )( *rng ) ];

#if VERBOSE >= 1
  cout << "[?: ";
  for ( auto it = k_nb.begin(); it != k_nb.end() - 1; ++it ) {
    cout << *it << " ";
  }
  cout << *( k_nb.end() - 1 )  << "] ";
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
  // assert( lat->chk_nn( hop.l, hop.k_pos ) );

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
