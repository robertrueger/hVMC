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

#include <vector>
#include <random>
#include <memory>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "macros.h"
#include "fptype.hpp"
#include "detwf.hpp"
#include "fpctrl.hpp"
#include "lattice.hpp"
#include "jastrow.hpp"
#include "econf.hpp"


class HubbardModelVMC
{
  protected:

    // ----- independent objects -----

    // random number generator
    const std::shared_ptr<std::mt19937> rng;

    // the underlying lattice
    const std::shared_ptr<Lattice> lat;

    // wavefunction and Jastrow
    const SingleParticleOrbitals detwf;
    const Jastrow v;

    // Hubbard model parameters
    const unsigned int update_hop_maxdist;
    const std::vector<fptype> t;
    const fptype U;


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

#ifdef USE_CBLAS
    // a temporary vectors that are large enough to hold one row/col of Wbu
    Eigen::VectorXfp tempWcol;
    Eigen::VectorXfp tempWrow;
#endif

    // buffer vector for X nearest neighbors
    // (in order to avoid allocating new ones all the time)
    std::vector<unsigned int> k_pos_Xnn;

    // floating point precision control
    const unsigned int updates_until_W_recalc, updates_until_T_recalc;
    unsigned int updates_since_W_recalc, updates_since_T_recalc;
    FPDevStat W_devstat, T_devstat;


    // ----- internal helper functions -----

    // function that performs a single Metropolis update
    bool metstep();

    // wrapper functions that updates/recalculates W/T after a successful hop
    void perform_W_update( const ElectronHop& hop );
    void perform_T_update( const ElectronHop& hop );

    // update and recalc functions for the internal objects
    void calc_new_W();
    void calc_qupdated_W( const ElectronHop& hop );
    Eigen::VectorXfp calc_new_T() const;
    Eigen::VectorXfp calc_qupdated_T( const ElectronHop& hop ) const;

    // functions to calculate the matrix D
    Eigen::MatrixXfp calc_Db() const;
    Eigen::MatrixXfp calc_Du() const;
    Eigen::MatrixXfp calc_Dd() const;


  public:

    HubbardModelVMC(
      const std::shared_ptr<std::mt19937>& rng_init,
      const std::shared_ptr<Lattice>& lat_init,
      const SingleParticleOrbitals& detwf_init,
      const Jastrow& v_init,
      unsigned int N_init,
      unsigned int update_hop_maxdist_init,
      const std::vector<fptype>& t_init, fptype U_init,
      fptype W_deviation_target_init,
      unsigned int updates_until_W_recalc_init,
      fptype T_deviation_target_init,
      unsigned int updates_until_T_recalc_init
    );

    // Monte Carlo step
    void mcs();

    // observable measurements
    fptype E_l();

    // floating point precision control
    FPDevStat get_W_devstat() const;
    FPDevStat get_T_devstat() const;

};

#endif // HUBBARD_MODEL_VMC_H_INCLUDED
