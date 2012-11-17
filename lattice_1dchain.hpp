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

#ifndef LATTICE_1DCHAIN_H_INCLUDED
#define LATTICE_1DCHAIN_H_INCLUDED

#include <vector>
#include <set>
#include <utility>
#include <algorithm>

#if VERBOSE >= 1
# include <iostream>
#endif

#include "macros.h"
#include "lattice.hpp"

class Lattice1DChain final : public Lattice {

  public:

    Lattice1DChain( unsigned int L_init );

    std::vector<unsigned int> get_Xnn( unsigned int l, unsigned int X ) const;

    IrreducibleIdxRel reduce_idxrel( unsigned int i, unsigned int j ) const;
    std::set<IrreducibleIdxRel> irreducible_idxrel_list() const;

};

#endif // LATTICE_1DCHAIN_H_INCLUDED
