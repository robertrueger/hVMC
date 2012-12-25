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

#ifndef ELECTRON_CONFIGURATION_H_INCLUDED
#define ELECTRON_CONFIGURATION_H_INCLUDED

#if VERBOSE >= 1 
# include <iostream>
#endif

#include <vector>
#include <random>
#include <stdexcept>

#include <CL/cl_platform.h>
#include <eigen3/Eigen/Eigen>

#include "macros.h"
#include "lattice.hpp"

enum ElectronOccupation_t {
  ELECTRON_OCCUPATION_EMPTY = 0,
  ELECTRON_OCCUPATION_FULL = 1
};


struct ElectronHop final {

  // id of the hopping electron
  const cl_uint k;

  // site that it hops to
  const cl_uint l;

  // position of electron k before the hop
  const cl_uint k_pos;

  // hop possible = site l unoccupied?
  const bool possible;

  ElectronHop( cl_uint k_init, cl_uint l_init,
               cl_uint k_pos_init, bool possible_init )
    : k( k_init ), l( l_init ),
      k_pos( k_pos_init ), possible( possible_init ) { }
};


class ElectronConfiguration final
{

  private:

    Lattice* const lat;
    const cl_uint electron_number;
    Eigen::Matrix<cl_uint, Eigen::Dynamic, 1> site_occ;
    std::vector<cl_uint> electron_pos;

    std::mt19937* const rng;

    // buffer vectors for nearest-neighbors
    // (in order to avoid allocating new ones all the time)
    std::vector<cl_uint> k_1nb, k_2nb, k_3nb;

    void reconstr_electron_pos();

  public:

    ElectronConfiguration(
      Lattice* const lat_init, cl_uint N_init,
      std::mt19937* rng_init
    );

    void init_from_raw_elpos( const std::vector<cl_uint>& raw_electron_pos );

    void distribute_random();

    ElectronHop propose_random_hop( cl_uint update_hop_maxdist );
    void do_hop( const ElectronHop& hop );

    cl_uint get_electron_pos( cl_uint k ) const;
    cl_uint get_site_occ( cl_uint l ) const;
    cl_uint N() const;
    cl_uint get_num_dblocc() const;

    std::vector<cl_uint> get_site_occ_raw() const;
    std::vector<cl_uint> get_electron_pos_raw() const;

};

#endif // ELECTRON_CONFIGURATION_H_INCLUDED
