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
#include "detwf.hpp"
#include "fpctrl.hpp"
#include "lattice.hpp"
#include "jastrow.hpp"
#include "pconf.hpp"
#include "wmatrix.hpp"
#include "tvector.hpp"


class HubbardModelVMC
{
  protected:

    // ----- independent objects -----

    // random number generator
    std::mt19937 rng;

    // the underlying lattice
    const std::shared_ptr<Lattice> lat;

    // wavefunction and Jastrow
    const DeterminantalWavefunction detwf;
    const Jastrow v;

    // Hubbard model parameters
    const unsigned int update_hop_maxdist;
    const std::vector<double> t;
    const double U;


    // ----- dependent and internal objects -----

    ParticleConfiguration pconf;

    WMatrix W;
    TVector T;

    // buffer vector for X nearest neighbors
    // (in order to avoid allocating new ones all the time)
    mutable std::vector<unsigned int> k_pos_Xnn;


    // ----- internal helper functions -----

    // function that performs a single Metropolis update
    bool metstep();


  public:

    HubbardModelVMC(
      const std::mt19937& rng_init,
      const std::shared_ptr<Lattice>& lat_init,
      const DeterminantalWavefunction& detwf_init,
      const Jastrow& v_init,
      unsigned int Ne_init,
      unsigned int update_hop_maxdist_init,
      const std::vector<double>& t_init, double U_init,
      double W_deviation_target,
      unsigned int updates_until_W_recalc,
      double T_deviation_target,
      unsigned int updates_until_T_recalc
    );

    // Monte Carlo step
    void mcs();

    // observable measurements
    double E_l() const;
    Eigen::VectorXd Delta_k( unsigned int optimizers ) const;
    double dblocc_dens() const;
    Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> n() const;
    Eigen::VectorXi s() const;

    // floating point precision control
    FPDevStat get_W_devstat() const;
    FPDevStat get_T_devstat() const;

};

#endif // HUBBARD_MODEL_VMC_H_INCLUDED
