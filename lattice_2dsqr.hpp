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

#ifndef LATTICE_2DSQUARE_H_INCLUDED
#define LATTICE_2DSQUARE_H_INCLUDED

#include <vector>
#include <set>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "macros.h"
#include "lattice.hpp"

class Lattice2DSquare : public Lattice {

  private:

    // side length of the square lattice of L sites
    unsigned int S;

    void get_1nn( unsigned int l, std::vector<unsigned int>* outbuf ) const;
    void get_2nn( unsigned int l, std::vector<unsigned int>* outbuf ) const;
    void get_3nn( unsigned int l, std::vector<unsigned int>* outbuf ) const;

    // functions to calculate positions on the lattice from spindices
    int x( Lattice::spindex l ) const;
    int y( Lattice::spindex l ) const;

    // functions to calculate position differences wrapped around the PBC
    int d( int p1, int p2 ) const;

  public:

    Lattice2DSquare( unsigned int L_init );

    void get_Xnn(
      spindex l, unsigned int X, std::vector<spindex>* outbuf
    ) const;

    irrspidxrel reduce_spidxrel( spindex i, spindex j ) const;
    std::set<irrspidxrel> get_all_irrspidxrels() const;
    irrspidxrel get_maxdist_irrspidxrel() const;

    Eigen::VectorXd r( index i, index j ) const;
    std::vector<Eigen::VectorXd> get_qvectors() const;

    double pairsym_modifier( optpairsym_t sym, spindex i, spindex j ) const;

    bool is_bipartite() const { return true; }
    unsigned int get_index_sublattice( index i ) const;
};

#endif // LATTICE_2DSQUARE_H_INCLUDED
