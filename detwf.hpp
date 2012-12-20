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

#ifndef DETERMINANTAL_WAVEFUNCTIONS_H_INCLUDED
#define DETERMINANTAL_WAVEFUNCTIONS_H_INCLUDED

# include <iostream>
#include <vector>

#include <CL/cl_platform.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>

#include "macros.h"
#include "fptype.hpp"
#include "lattice.hpp"


struct SingleParticleOrbitals {

  // the actual orbitals
  const Eigen::MatrixXfp orbitals;

  // spin symmetric and orbitals spin eigenstates
  const bool ssym;

  SingleParticleOrbitals(
    const Eigen::MatrixXfp& orbitals_init, bool ssym_init
  ) : orbitals( orbitals_init ), ssym( ssym_init ) { }
};


SingleParticleOrbitals wf_tight_binding(
    const std::vector<cl_fptype>& t,
    cl_uint N, Lattice* lat
);

bool check_openshell( const Eigen::VectorXfp& E, cl_uint N );

#endif // DETERMINANTAL_WAVEFUNCTIONS_H_INCLUDED
