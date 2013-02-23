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
#include "fptype.hpp"
#include "lattice.hpp"
#include "varparam.hpp"


class Jastrow
{

  private:

    const std::shared_ptr<Lattice> lat;

    std::vector<fptype> idxrel_expv;

  public:

    Jastrow(
      const std::shared_ptr<Lattice>& lat_init,
      const Eigen::VectorXfp& v_init
    );

    fptype operator()( unsigned int i, unsigned int j ) const;
    fptype exp( unsigned int i, unsigned int j ) const;
    fptype exp_onsite() const;
    void set( unsigned int i, unsigned int j, fptype v_new  );

};

#endif // JASTROW_H_INCLUDED
