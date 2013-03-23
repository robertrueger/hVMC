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

#ifndef WMATRIX_H_INCLUDED
#define WMATRIX_H_INCLUDED

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "fptype.hpp"
#include "econf.hpp"
#include "lattice.hpp"
#include "detwf.hpp"
#include "fpctrl.hpp"


class WMatrix final
{
  protected:

    const Lattice* const lat;
    const SingleParticleOrbitals& detwf;
    const ElectronConfiguration& econf;

    Eigen::MatrixXfp  Wbu_1;
    Eigen::MatrixXfp  Wbu_2;
    Eigen::MatrixXfp* Wbu_active;
    Eigen::MatrixXfp* Wbu_inactive;

    Eigen::MatrixXfp  Wd_1;
    Eigen::MatrixXfp  Wd_2;
    Eigen::MatrixXfp* Wd_active;
    Eigen::MatrixXfp* Wd_inactive;

#ifdef USE_CBLAS
    // a temporary vectors that are large enough to hold one row/col of Wbu
    Eigen::VectorXfp tempWcol;
    Eigen::VectorXfp tempWrow;
#endif

    const unsigned int updates_until_recalc;
    unsigned int updates_since_recalc;
    FPDevStat devstat;

    // functions to calculate the matrix D
    Eigen::MatrixXfp calc_Db() const;
    Eigen::MatrixXfp calc_Du() const;
    Eigen::MatrixXfp calc_Dd() const;

    void calc_new();
    void calc_qupdated( const ElectronHop& hop );

  public:

    WMatrix(
      const Lattice* lat_init,
      const SingleParticleOrbitals& detwf_init,
      const ElectronConfiguration& econf_init,
      fptype deviation_target,
      unsigned int updates_until_recalc_init
    );

    bool init_and_check();
    void update( const ElectronHop& hop );

    fptype operator()( unsigned int i, unsigned int j ) const;

    FPDevStat get_devstat() const;

};

#endif // WMATRIX_H_INCLUDED
