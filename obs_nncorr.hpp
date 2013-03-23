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

#ifndef OBS_DENSITY_DENSITY_CORRELATION_H_INCLUDED
#define OBS_DENSITY_DENSITY_CORRELATION_H_INCLUDED

#include "obs.hpp"

#include <vector>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "fptype.hpp"


class ObservableDensityDensityCorrelation final : public Observable
{
  private:

    Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> nncorr_sum;
    unsigned int this_bin_num_measurements;

    std::vector<Eigen::MatrixXfp> nncorr_binmeans;

  public:

    ObservableDensityDensityCorrelation();

    void measure( HubbardModelVMC& model, ObservableCache& cache );

    void completebin();

    void collect_and_write_results(
      const boost::mpi::communicator& mpicomm,
      MCCResults& results
    ) const;

    void send_results_to_master( const boost::mpi::communicator& mpicomm ) const;
};

#endif // OBS_DENSITY_DENSITY_CORRELATION_H_INCLUDED