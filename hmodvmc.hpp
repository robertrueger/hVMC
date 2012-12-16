/*
 * Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
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

#include <CL/cl_platform.h>

#include "fptype.hpp"
#include "fpctrl.hpp"


class HubbardModelVMC
{
  public:

    virtual ~HubbardModelVMC() { }

    // Monte Carlo step
    virtual void mcs() = 0;
    virtual void equilibrate( cl_uint N_mcs_equil ) = 0;

    // observable measurements
    virtual cl_fptype E_l() = 0;
    virtual cl_ulong mctime() const = 0;

    // floating point precision control
    virtual FPDevStat get_W_devstat() const = 0;
    virtual FPDevStat get_T_devstat() const = 0;

};

#endif // HUBBARD_MODEL_VMC_H_INCLUDED
