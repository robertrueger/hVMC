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

#ifndef OBS_CORRELATION_H_INCLUDED
#define OBS_CORRELATION_H_INCLUDED

#include "obs.hpp"

#include <boost/mpi/communicator.hpp>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "modman.hpp"
#include "mccresults.hpp"

class ObservableCorrelation : public Observable
{
  protected:

    Eigen::MatrixXd thisbin_sum;
    unsigned int thisbin_count;

    Eigen::MatrixXd binmean_sum;
    unsigned int binmean_count;

    // methods that are specialized to the actual correlation being measured
    virtual Eigen::VectorXd get_current(
      const ModelManager& model, ObservableCache& cache
    ) const = 0;
    virtual void save_to_results(
      const Eigen::MatrixXd& corrresult, MCCResults& results
    ) const = 0;

  public:

    ObservableCorrelation( unsigned int L );

    void measure( const ModelManager& model, ObservableCache& cache );

    void completebin();

    void collect_and_write_results(
      const boost::mpi::communicator& mpicomm,
      MCCResults& results
    ) const;

    void send_results_to_master( const boost::mpi::communicator& mpicomm ) const;
};

#endif // OBS_CORRELATION_H_INCLUDED
