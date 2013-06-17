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

#ifndef LATTICE_2DSQUARE_2LAYER_H_INCLUDED
#define LATTICE_2DSQUARE_2LAYER_H_INCLUDED

#include <vector>
#include <set>
#include <random>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "macros.h"
#include "lattice.hpp"


class Lattice2DSquare2Layer : public Lattice {

  private:

    // number of sites per layer
    const unsigned int L_layer;

    // side length of the square lattice of L sites
    const unsigned int S;

    // first and second nearest neighbor in the same plane
    void get_1nn( unsigned int l, std::vector<unsigned int>* outbuf ) const;
    void get_2nn( unsigned int l, std::vector<unsigned int>* outbuf ) const;
    // corresponding site in the other plane (treated as third n.n.)
    void get_opn( unsigned int l, std::vector<unsigned int>* outbuf ) const;

    // functions to calculate positions on the lattice from spindices
    int x( Lattice::spindex l ) const;
    int y( Lattice::spindex l ) const;
    int z( Lattice::spindex l ) const;

    // functions to calculate position differences wrapped around the PBC
    int d( int p1, int p2 ) const;

  public:

    Lattice2DSquare2Layer( unsigned int L_init );

    void get_Xnn(
      spindex l, unsigned int X, std::vector<spindex>* outbuf
    ) const;

    irridxrel reduce_idxrel( spindex i, spindex j ) const;
    std::set<irridxrel> get_all_irridxrels() const;
    irridxrel get_maxdist_irridxrel() const;

    Eigen::VectorXd r( index i, index j ) const;
    bool include_r_in_ssfac( index i, index j ) const;
    std::vector<Eigen::VectorXd> get_qvectors() const;

    double pairsym_modifier( optpairsym_t sym, spindex i, spindex j ) const;

    Eigen::VectorXi get_random_spindex_occ(
      unsigned int Npu, unsigned int Npd, std::mt19937& rng ) const;
};

#endif // LATTICE_2DSQUARE_2LAYER_H_INCLUDED
