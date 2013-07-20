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

#ifndef JASTROW_H_INCLUDED
#define JASTROW_H_INCLUDED

#include <vector>
#include <memory>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "macros.h"
#include "lattice.hpp"
#include "varparam.hpp"


class Jastrow
{

  private:

    const std::shared_ptr<Lattice> lat;

    std::vector<double> iir_v;

    unsigned int num_vpar;
    std::vector<unsigned int> iir_vparnum;

  public:

    Jastrow(
      const std::shared_ptr<Lattice>& lat_init,
      const Eigen::VectorXd& v_init
    );

    double operator()( Lattice::spindex i, Lattice::spindex j ) const;
    double onsite() const;
    void set( Lattice::spindex i, Lattice::spindex j, double v_new  );

    unsigned int get_num_vpar() const;
    unsigned int get_vparnum( Lattice::irridxrel iir ) const;
};

#endif // JASTROW_H_INCLUDED
