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

#ifndef HUBBARD_MODEL_VMC_CPU_H_INCLUDED
#define HUBBARD_MODEL_VMC_CPU_H_INCLUDED

#if VERBOSE >= 1
# include <iostream>
#endif

#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

#include <CL/cl_platform.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

#include "macros.h"
#include "fptype.hpp"
#include "detwf.hpp"
#include "fpctrl.hpp"
#include "lattice.hpp"
#include "jastrow.hpp"
#include "econf.hpp"
#include "hmodvmc.hpp"


class HubbardModelVMC_CPU final : public HubbardModelVMC
{
  private:

    // forbid copy construction and assignment
    HubbardModelVMC_CPU( const HubbardModelVMC_CPU& other );
    HubbardModelVMC_CPU& operator=( const HubbardModelVMC_CPU& other );
    // (both are not implemented)

  protected:

    // ----- independent objects -----

    // random number generator
    std::mt19937 rng;

    // the underlying lattice
    Lattice* const lat;

    // wavefunction and Jastrow
    const SingleParticleOrbitals M;
    const Jastrow v;

    // Hubbard model parameters
    const cl_uint update_hop_maxdist;
    const std::vector<cl_fptype> t;
    const cl_fptype U;


    // ----- dependent and internal objects -----

    ElectronConfiguration econf;

    Eigen::MatrixXfp  Wbu_1;
    Eigen::MatrixXfp  Wbu_2;
    Eigen::MatrixXfp* Wbu_active;
    Eigen::MatrixXfp* Wbu_inactive;

    Eigen::MatrixXfp  Wd_1;
    Eigen::MatrixXfp  Wd_2;
    Eigen::MatrixXfp* Wd_active;
    Eigen::MatrixXfp* Wd_inactive;

    Eigen::VectorXfp T;

    cl_ulong completed_mcsteps;

    // floating point precision control
    const cl_uint updates_until_W_recalc, updates_until_T_recalc;
    cl_uint updates_since_W_recalc, updates_since_T_recalc;
    FPDevStat W_devstat, T_devstat;


    // ----- internal helper functions -----

    // function that performs a single Metropolis update
    bool metstep();

    // wrapper functions that updates/recalculates W/T after a successful hop
    void perform_W_update( const ElectronHop& hop );
    void perform_T_update( const ElectronHop& hop );

    // update and recalc functions for the internal objects
    void calc_new_W();
    void calc_qupdated_Wbu( const ElectronHop& hop );
    void calc_qupdated_Wd(  const ElectronHop& hop );
    Eigen::VectorXfp calc_new_T() const;
    Eigen::VectorXfp calc_qupdated_T( const ElectronHop& hop ) const;

    // functions to calculate the matrix D
    Eigen::MatrixXfp calc_Db() const;
    Eigen::MatrixXfp calc_Du() const;
    Eigen::MatrixXfp calc_Dd() const;


  public:

    HubbardModelVMC_CPU(
      std::mt19937 rng_init,
      Lattice* const lat_init,
      const SingleParticleOrbitals& M_init,
      const Jastrow& v_init,
      cl_uint N_init,
      cl_uint update_hop_maxdist_init,
      const std::vector<cl_fptype>& t_init, cl_fptype U_init,
      cl_fptype W_deviation_target_init,
      cl_uint updates_until_W_recalc_init,
      cl_fptype T_deviation_target_init,
      cl_uint updates_until_T_recalc_init
    );

    ~HubbardModelVMC_CPU();

    // Monte Carlo step
    void mcs();
    void equilibrate( cl_uint N_mcs_equil );

    // observable measurements
    cl_fptype E_l() const;
    cl_ulong mctime() const;

    // floating point precision control
    FPDevStat get_W_devstat() const;
    FPDevStat get_T_devstat() const;

};

#endif // HUBBARD_MODEL_VMC_CPU_H_INCLUDED
