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

#ifndef UTILS_CL_H_INCLUDED
#define UTILS_CL_H_INCLUDED

#include <iostream>
#include <vector>
#include <string>
#include <iterator>
#include <stdexcept>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>


void print_clinfo();

cl::Context clcontext_setup( cl_uint pl_id, cl_uint dev_id );

#endif // UTILS_CL_H_INCLUDED
