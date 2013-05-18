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

#ifndef OBS_SPIN_SPIN_CORRELATION_H_INCLUDED
#define OBS_SPIN_SPIN_CORRELATION_H_INCLUDED

#include "obs.hpp"
#include "obs_corr.hpp"

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "hmodvmc.hpp"
#include "mccresults.hpp"


class ObservableSpinSpinCorrelation final : public ObservableCorrelation
{
  protected:

    Eigen::MatrixXd get_current(
      const HubbardModelVMC& model, ObservableCache& cache
    ) const;

    void save_to_results(
      const Eigen::MatrixXd& corrresult, MCCResults& results
    ) const;

  public:

    ObservableSpinSpinCorrelation( unsigned int L )
      : ObservableCorrelation( L ) { };
};

#endif // OBS_SPIN_SPIN_CORRELATION_H_INCLUDED
