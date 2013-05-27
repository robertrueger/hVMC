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

#ifndef OBSERVABLES_H_INCLUDED
#define OBSERVABLES_H_INCLUDED

#include <boost/mpi/communicator.hpp>
#include <boost/optional.hpp>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "hmodvmc.hpp"
#include "mccresults.hpp"


enum observables_t {
  OBSERVABLE_E,
  OBSERVABLE_DELTAK,
  OBSERVABLE_DELTAK_DELTAKPRIME,
  OBSERVABLE_DELTAK_E,
  OBSERVABLE_DOUBLE_OCCUPANCY_DENSITY,
  OBSERVABLE_DENSITY_DENSITY_CORRELATION,
  OBSERVABLE_SPIN_SPIN_CORRELATION,
  OBSERVABLE_PARTICLE_CONFIGURATIONS
};


struct ObservableCache
{
  boost::optional<double> E;
  boost::optional<Eigen::VectorXd> DeltaK;
  boost::optional<double> dblocc;
  boost::optional< Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> > n;
  boost::optional< Eigen::VectorXi > s;

  void clear() {
    E = boost::none;
    DeltaK = boost::none;
    dblocc = boost::none;
    n = boost::none;
    s = boost::none;
  }
};


class Observable
{
  public:

    virtual ~Observable() { };

    virtual void measure(
      const HubbardModelVMC& model,
      ObservableCache& cache
    ) = 0;

    virtual void completebin() = 0;

    virtual void collect_and_write_results(
      const boost::mpi::communicator& mpicomm,
      MCCResults& results
    ) const = 0;
    virtual void send_results_to_master(
      const boost::mpi::communicator& mpicomm
    ) const = 0;
};

#endif // OBSERVABLES_H_INCLUDED
