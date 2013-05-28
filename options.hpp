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

#ifndef OPTIONS_H_INCLUDED
#define OPTIONS_H_INCLUDED

#include <boost/program_options/variables_map.hpp>

enum optmode_t {
  OPTION_MODE_OPTIMIZATION,
  OPTION_MODE_SIMULATION,
  OPTION_MODE_ANALYSIS
};

enum optpairsym_t {
  OPTION_PAIRING_SYMMETRY_SWAVE,
  OPTION_PAIRING_SYMMETRY_DWAVE,
  OPTION_PAIRING_SYMMETRY_DWAVE_TWISTED
};

typedef boost::program_options::variables_map Options;

Options read_options( int argc, char* argv[], bool is_master );

#endif // OPTIONS_H_INCLUDED
