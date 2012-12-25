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

#ifndef HUBBARD_MODEL_VMC_CL_H_INCLUDED
#define HUBBARD_MODEL_VMC_CL_H_INCLUDED

#if VERBOSE >= 1
# include <iostream>
#endif

#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

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
#include "utils.hpp"


class HubbardModelVMC_CL final : public HubbardModelVMC
{
  private:

    // forbid copy construction and assignment
    HubbardModelVMC_CL( const HubbardModelVMC_CL& other );
    HubbardModelVMC_CL& operator=( const HubbardModelVMC_CL& other );
    // (both are not implemented)

  protected:

    // ----- OpenCL objects -----

    cl::Context clCtx;
    cl::Program clPrg;
    cl::CommandQueue clQ;


    // ----- independent objects -----

    // random number generator
    std::mt19937 rng;

    // the underlying lattice
    Lattice* const lat;

    // wavefunction and Jastrow
    const SingleParticleOrbitals M;
    const Jastrow v;

    // Hubbard model parameters
    const std::vector<cl_fptype> t;
    const cl_fptype U;


    // ----- dependent and internal objects ----- 
    
    // device buffers for the jastrow factors and the hopping parameters + U
    cl::Buffer devexpv;
    cl::Buffer devUt;

    // device buffer for the current electron hop
    cl::Buffer devehop;

    // electron configuration on the host
    ElectronConfiguration econf;

    // electron configuration device buffers
    cl::Buffer deveconf_site_occ;
    cl::Buffer deveconf_electron_pos;

    // host buffers for the matrix W
    Eigen::MatrixXfp Wbu;
    Eigen::MatrixXfp Wd;
    Eigen::MatrixXfp Wbu_fromdev;
    Eigen::MatrixXfp Wd_fromdev;

    // device buffers for matrix W
    cl::Buffer  devWbu_1;
    cl::Buffer  devWbu_2;
    cl::Buffer* devWbu_active;
    cl::Buffer* devWbu_inactive;

    cl::Buffer  devWd_1;
    cl::Buffer  devWd_2;
    cl::Buffer* devWd_active;
    cl::Buffer* devWd_inactive;

    // host buffers for vector T
    Eigen::VectorXfp T;
    Eigen::VectorXfp T_fromdev;

    // device buffers for vector T
    cl::Buffer  devT_1;
    cl::Buffer  devT_2;
    cl::Buffer* devT_active;
    cl::Buffer* devT_inactive;

    // buffer for each electron's part of the local energy
    Eigen::VectorXfp E_l_elbuf;
    cl::Buffer devE_l_elbuf;

    cl_ulong completed_mcsteps;

    // floating point precision control
    const cl_uint updates_until_W_recalc, updates_until_T_recalc;
    cl_uint updates_since_W_recalc, updates_since_T_recalc;
    FPDevStat W_devstat, T_devstat;


    // ----- internal helper functions -----
    
    // sync electronic configuration on the host with the device
    void sync_econf_down();

    // function that performs a single Metropolis update
    void metstep();
    cl::Kernel clK_hop;

    // wrapper functions that update/recalculate W/T after a successful hop
    void perform_W_update();
    void perform_T_update();

    // update and recalc functions for the internal objects
    void calc_new_W();
    void calc_qupdated_W();
    cl::Kernel clK_update_devWbu;
    cl::Kernel clK_update_devWd;

    void calc_new_T();
    void calc_qupdated_T();
    cl::Kernel clK_update_devT;

    // local energy calculation kernel
    cl::Kernel clK_calc_E_l;

    // functions to calculate the matrix D
    Eigen::MatrixXfp calc_Db() const;
    Eigen::MatrixXfp calc_Du() const;
    Eigen::MatrixXfp calc_Dd() const;


  public:

    HubbardModelVMC_CL(
      const cl::Context& clCtx_init,
      std::mt19937 rng_init,
      Lattice* const lat_init,
      const SingleParticleOrbitals& M_init,
      const Jastrow& v_init,
      cl_uint N_init,
      const std::vector<cl_fptype>& t_init, cl_fptype U_init,
      cl_fptype W_deviation_target_init,
      cl_uint updates_until_W_recalc_init,
      cl_fptype T_deviation_target_init,
      cl_uint updates_until_T_recalc_init
    );

    ~HubbardModelVMC_CL();

    // Monte Carlo step
    void mcs();
    void equilibrate( cl_uint N_mcs_equil );

    // observable measurements
    cl_fptype E_l();
    cl_ulong mctime() const;

    // floating point precision control
    FPDevStat get_W_devstat() const;
    FPDevStat get_T_devstat() const;

};

#endif // HUBBARD_MODEL_VMC_CL_H_INCLUDED
