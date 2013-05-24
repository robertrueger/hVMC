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

#ifndef PARTICLE_CONFIGURATION_H_INCLUDED
#define PARTICLE_CONFIGURATION_H_INCLUDED

#include <vector>
#include <random>
#include <memory>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "macros.h"
#include "lattice.hpp"


enum ParticleOccupation_t {
  PARTICLE_OCCUPATION_EMPTY = 0,
  PARTICLE_OCCUPATION_FULL = 1
};


struct ParticleHop final {

  // id of the hopping particle
  const unsigned int k;

  // site that it hops to
  const Lattice::spindex l;

  // position of particle k before the hop
  const Lattice::spindex k_pos;

  // hop possible = site l unoccupied?
  const bool possible;

  ParticleHop( Lattice::spindex k_init, Lattice::spindex l_init,
               Lattice::spindex k_pos_init, bool possible_init )
    : k( k_init ), l( l_init ),
      k_pos( k_pos_init ), possible( possible_init ) { }
};


class ParticleConfiguration final
{

  private:

    const std::shared_ptr<Lattice> lat;

  public:

    // BEWARE: particle-hole-transformation!!
    // -> spin up particles are spin up electrons
    // -> spin down particles are spin down electron holes
    // total number of ...
    const unsigned int Ne; // ... electrons
    const unsigned int Npu; // ... spin up particles
    const unsigned int Npd; /// ... spin down particles
    const unsigned int Np; // ... particles
    Eigen::VectorXi spindex_occ;
    std::vector<Lattice::spindex> particlenum_pos;

  private:

    std::mt19937& rng;

    // buffer vectors for nearest-neighbors
    // (in order to avoid allocating new ones all the time)
    mutable std::vector<Lattice::spindex> k_1nb, k_2nb, k_3nb;

    void reconstr_particlenum_pos();

  public:

    ParticleConfiguration(
      const std::shared_ptr<Lattice>& lat_init, unsigned int Ne_init,
      std::mt19937& rng_init
    );

    void distribute_random();

    ParticleHop propose_random_hop( unsigned int update_hop_maxdist ) const;
    void do_hop( const ParticleHop& hop );

    const Eigen::VectorXi& get_spindex_occ() const;
    const std::vector<Lattice::spindex>& get_particlenum_pos() const;
};


#endif // PARTICLE_CONFIGURATION_H_INCLUDED
