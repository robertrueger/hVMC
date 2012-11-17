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

#ifndef HUBBARD_MODEL_VMC_H_INCLUDED
#define HUBBARD_MODEL_VMC_H_INCLUDED

#if VERBOSE >= 1
# include <iostream>
#endif

#include <vector>
#include <random>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

#include "macros.h"
#include "fptype.hpp"
#include "lattice.hpp"
#include "jastrow.hpp"
#include "econf.hpp"


class HubbardModelVMC final
{
  private:

    // forbid copy construction and assignment
    HubbardModelVMC( const HubbardModelVMC& other );
    HubbardModelVMC& operator=( const HubbardModelVMC& other );
    // (both are not implemented)

  protected:

    // ----- independent objects -----

    // random number generator
    std::mt19937 rng;

    // the underlying lattice
    Lattice* const lat;

    // wavefunction and Jastrow
    const Eigen::MatrixXfp M;
    const Jastrow v;

    // Hubbard model parameters
    const unsigned int update_hop_maxdist;
    const std::vector<fptype> t;
    const fptype U;


    // ----- dependent and internal objects -----

    // Monte Carlo cycle
    ElectronConfiguration econf;
    Eigen::MatrixXfp W;
    Eigen::VectorXfp T;
    unsigned long int completed_mcsteps;

    // number of updates until recalc of W and T
    const unsigned int updates_until_WT_recalc;
    unsigned int updates_since_WT_recalc;


    // ----- internal helper functions -----

    // function that performs a single Metropolis update
    bool metstep();

    // wrapper function that updates/recalculates WT after a successful hop
    void perform_WT_update( const ElectronHop& hop );

  public:

    HubbardModelVMC(
      std::mt19937 rng_init,
      Lattice* const lat_init,
      const Eigen::MatrixXfp M_init,
      const Jastrow& v_init,
      unsigned int N_init,
      unsigned int update_hop_maxdist_init,
      const std::vector<fptype>& t_init, fptype U_init,
      unsigned int updates_until_WT_recalc_init
    );

    ~HubbardModelVMC();

    // Monte Carlo step
    void mcs();
    void equilibrate( unsigned int N_mcs_equil );

    // helper functions (can be public, because they are const)
    Eigen::MatrixXfp calc_D() const;
    Eigen::MatrixXfp calc_new_W() const;
    Eigen::MatrixXfp calc_updated_W( const ElectronHop& hop ) const;
    Eigen::VectorXfp calc_new_T() const;
    Eigen::VectorXfp calc_updated_T( const ElectronHop& hop ) const;

    // observable measurements
    fptype E_l() const;
    unsigned long int mctime() const;

};

#endif // HUBBARD_MODEL_VMC_H_INCLUDED
