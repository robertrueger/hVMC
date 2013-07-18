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

#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED

#include <vector>
#include <set>
#include <random>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "options.hpp"
#include "macros.h"


class Lattice {

  public:

    // the lattice type
    enum class type {
        chain1d,
        square2d,
        square2d2layer
    };

    // the number of lattice sites
    const unsigned int L;

    Lattice( unsigned int L_init ) : L( L_init ) { }
    virtual ~Lattice() { }

    // typedefs for the different types of indices
    typedef unsigned int index;
    typedef unsigned int spindex;
    enum class spindex_type { up, down };
    typedef unsigned int irridxrel;

    // convenience functions for the indices
    index get_index_from_spindex( spindex l ) const;
    spindex get_linked_spindex( spindex l ) const;
    spindex_type get_spindex_type( spindex l ) const;

    // method to find the neighboring sites
    std::vector<spindex> get_Xnn( spindex l, unsigned int X ) const;
    virtual void get_Xnn(
      spindex l, unsigned int X, std::vector<spindex>* outbuf
    ) const = 0;

    // everything related to reduced index relations
    virtual irridxrel reduce_idxrel( spindex i, spindex j ) const = 0;
    virtual std::set<irridxrel> get_all_irridxrels() const = 0;
    virtual irridxrel get_maxdist_irridxrel() const = 0;

    // lattice geometry and relevant reciprocal lattice vectors
    virtual Eigen::VectorXd r( index i, index j ) const = 0;
    virtual bool include_r_in_ssfac( index i, index j ) const;
    virtual std::vector<Eigen::VectorXd> get_qvectors() const = 0;

    // pairing symmetry modifying factor
    virtual double pairsym_modifier(
        optpairsym_t sym, spindex i, spindex j
    ) const = 0;

    // bipartite lattices and lattice site types
    virtual bool is_bipartite() const { return false; }
    virtual unsigned int get_index_sublattice( index ) const { return 0; }

    // method to generate a random particle distribution
    // (this is part of the lattice class because there might be special
    // requirements for a random particle distribution, e.g.: the two layered
    // Hubbard model requires each plane to have the same number of particles)
    virtual Eigen::VectorXi get_random_spindex_occ(
      unsigned int Npu, unsigned int Npd, std::mt19937& rng
    ) const;
};

#endif // LATTICE_H_INCLUDED
