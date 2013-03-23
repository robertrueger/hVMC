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
#include "wmatrix.hpp"


class HubbardModelVMC final
{
  protected:

    // ----- independent objects -----

    // random number generator
    std::mt19937 rng;

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

    WMatrix W;
    Eigen::VectorXfp T;

    // buffer vector for X nearest neighbors
    // (in order to avoid allocating new ones all the time)
    std::vector<unsigned int> k_pos_Xnn;

    // floating point precision control
    const unsigned int updates_until_T_recalc;
    unsigned int updates_since_T_recalc;
    FPDevStat T_devstat;


    // ----- internal helper functions -----

    // function that performs a single Metropolis update
    bool metstep();

    // wrapper functions that updates/recalculates W/T after a successful hop
    void perform_T_update( const ElectronHop& hop );

    // update and recalc functions for the internal objects
    Eigen::VectorXfp calc_new_T() const;
    Eigen::VectorXfp calc_qupdated_T( const ElectronHop& hop ) const;



  public:

    HubbardModelVMC(
      const std::mt19937& rng_init,
      const std::shared_ptr<Lattice>& lat_init,
      const SingleParticleOrbitals& detwf_init,
      const Jastrow& v_init,
      unsigned int N_init,
      unsigned int update_hop_maxdist_init,
      const std::vector<fptype>& t_init, fptype U_init,
      fptype W_deviation_target,
      unsigned int updates_until_W_recalc,
      fptype T_deviation_target_init,
      unsigned int updates_until_T_recalc_init
    );

    // Monte Carlo step
    void mcs();

    // observable measurements
    fptype E_l();
    Eigen::VectorXfp Delta_k() const;
    fptype dblocc_dens() const;
    Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> n() const;

    // floating point precision control
    FPDevStat get_W_devstat() const;
    FPDevStat get_T_devstat() const;

};

#endif // HUBBARD_MODEL_VMC_H_INCLUDED
