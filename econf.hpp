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

#ifndef ELECTRON_CONFIGURATION_H_INCLUDED
#define ELECTRON_CONFIGURATION_H_INCLUDED

#include <vector>
#include <random>
#include <memory>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "macros.h"
#include "lattice.hpp"


enum ElectronOccupation_t {
  ELECTRON_OCCUPATION_EMPTY = 0,
  ELECTRON_OCCUPATION_FULL = 1
};


struct ElectronHop final {

  // id of the hopping electron
  const unsigned int k;

  // site that it hops to
  const unsigned int l;

  // position of electron k before the hop
  const unsigned int k_pos;

  // hop possible = site l unoccupied?
  const bool possible;

  ElectronHop( unsigned int k_init, unsigned int l_init,
               unsigned int k_pos_init, bool possible_init )
    : k( k_init ), l( l_init ),
      k_pos( k_pos_init ), possible( possible_init ) { }
};


class ElectronConfiguration final
{

  private:

    const std::shared_ptr<Lattice> lat;
    const unsigned int electron_number;
    Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> site_occ;
    std::vector<unsigned int> electron_pos;

    const std::shared_ptr<std::mt19937> rng;

    // buffer vectors for nearest-neighbors
    // (in order to avoid allocating new ones all the time)
    std::vector<unsigned int> k_1nb, k_2nb, k_3nb;

    void reconstr_electron_pos();

  public:

    ElectronConfiguration(
      const std::shared_ptr<Lattice>& lat_init, unsigned int N_init,
      const std::shared_ptr<std::mt19937>& rng_init
    );

    void distribute_random();

    ElectronHop propose_random_hop( unsigned int update_hop_maxdist );
    void do_hop( const ElectronHop& hop );

    unsigned int get_electron_pos( unsigned int k ) const;
    unsigned int get_site_occ( unsigned int l ) const;
    unsigned int N() const;
    unsigned int get_num_dblocc() const;

};

#endif // ELECTRON_CONFIGURATION_H_INCLUDED
