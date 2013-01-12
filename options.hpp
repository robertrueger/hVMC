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

#ifndef OPTIONS_H_INCLUDED
#define OPTIONS_H_INCLUDED

#include <boost/program_options/variables_map.hpp>
#include <boost/mpi/communicator.hpp>

#include "macros.h"

typedef boost::program_options::variables_map Options;

Options read_options(
  int argc, char* argv[],
  const boost::mpi::communicator& mpiflock
);

#endif // OPTIONS_H_INCLUDED
