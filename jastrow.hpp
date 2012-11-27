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

#ifndef JASTROW_H_INCLUDED
#define JASTROW_H_INCLUDED

#if VERBOSE >= 1
# include <iostream>
#endif

#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

#include "macros.h"
#include "fptype.hpp"
#include "lattice.hpp"


class Jastrow final
{

  private:

    Lattice* const lat;

    std::vector< std::vector<float> > idxrel_expv;

  public:

    Jastrow( Lattice* lat_init );

    void randomize( fptype min, fptype max, std::mt19937* rng );

    fptype operator()( unsigned int i, unsigned int j ) const;
    fptype exp( unsigned int i, unsigned int j ) const;
    void set( unsigned int i, unsigned int j, fptype v_new  );

};

#endif // JASTROW_H_INCLUDED
