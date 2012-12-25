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
  cl_uint electron_number_init,
  mt19937* rng_init )
  : lat( lat_init ), electron_number( electron_number_init ),
    site_occ(
      Eigen::Matrix<cl_uint, Eigen::Dynamic, 1>::Zero( 2 * lat_init->L )
    ),
    electron_pos( std::vector<cl_uint>( electron_number_init ) ),
    rng( rng_init )
{
  distribute_random();
}



void ElectronConfiguration::init_from_raw_elpos(
  const std::vector<cl_uint>& raw_electron_pos )
{
  assert( raw_electron_pos.size() == electron_number );
  electron_pos = raw_electron_pos;
  for ( cl_uint i = 0; i < electron_pos.size(); ++i ) {
    site_occ( electron_pos[ i ] ) = ELECTRON_OCCUPATION_FULL;
  }
}



void ElectronConfiguration::reconstr_electron_pos()
{
  electron_pos.clear();
  for ( cl_uint l = 0; l < 2 * lat->L; ++l ) {
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
  site_occ = Eigen::Matrix<cl_uint, Eigen::Dynamic, 1>::Zero( 2 * lat->L );

  // randomly distribute L/2 electrons per spin direction
  while ( site_occ.head( lat->L ).sum() < electron_number / 2 ) {
    site_occ[ uniform_int_distribution<cl_uint>( 0, lat->L - 1 )( *rng ) ]
      = ELECTRON_OCCUPATION_FULL;
  }
  while ( site_occ.tail( lat->L ).sum() < electron_number / 2 ) {
    site_occ[ uniform_int_distribution<cl_uint>( 0, lat->L - 1 )( *rng )
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
  cl_uint update_hop_maxdist )
{
  // hop the kth electron
  const cl_uint k
    = uniform_int_distribution<cl_uint>( 0, electron_number - 1 )( *rng );

  // find the position of the xth electron
  const cl_uint k_pos = electron_pos[k];

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

  const cl_uint nb_number
    = uniform_int_distribution<cl_uint>
      ( 0,
        k_1nb.size() + k_2nb.size() + k_3nb.size() - 1
      )( *rng );

  cl_uint l;
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



cl_uint ElectronConfiguration::get_electron_pos( cl_uint k ) const
{
  return electron_pos[ k ];
}



cl_uint ElectronConfiguration::get_site_occ( cl_uint l ) const
{
  return site_occ[ l ];
}



cl_uint ElectronConfiguration::N() const
{
  return electron_number;
}



cl_uint ElectronConfiguration::get_num_dblocc() const
{
  return ( site_occ.head( lat->L ).array() * site_occ.tail( lat->L ).array() ).sum();
}



vector<cl_uint> ElectronConfiguration::get_site_occ_raw() const
{
  vector<cl_uint> rawdat( site_occ.size() );
  for ( cl_uint i = 0; i < site_occ.size(); ++i ) {
    rawdat.at( i ) = site_occ( i );
  }
  return rawdat;
}



vector<cl_uint> ElectronConfiguration::get_electron_pos_raw() const
{
  return electron_pos;
}
