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

#ifndef MCCRUN_H_INCLUDED
#define MCCRUN_H_INCLUDED

#include <set>

#include <boost/mpi/communicator.hpp>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "options.hpp"
#include "mccresults.hpp"
#include "obs.hpp"


MCCResults mccrun_master(
  const Options& opts, const Eigen::VectorXd& vpar, unsigned int num_bins,
  const std::set<observables_t>& obs, const boost::mpi::communicator& mpicomm,
  boost::optional<const Eigen::VectorXi&> pconf_init
    = boost::optional<const Eigen::VectorXi&>()
);

void mccrun_slave(
  const Options& opts, const Eigen::VectorXd& vpar,
  const std::set<observables_t>& obs, const boost::mpi::communicator& mpicomm,
  boost::optional<const Eigen::VectorXi&> pconf_init
    = boost::optional<const Eigen::VectorXi&>()
);

#endif // MCCRUN_H_INCLUDED
