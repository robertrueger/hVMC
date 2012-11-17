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

#ifndef SIMRESULTS_H_INCLUDED
#define SIMRESULTS_H_INCLUDED

#include "macros.h"
#include "fptype.hpp"


struct BasicSimResults final {

  bool success;

  fptype E;
  fptype sigma_E;

  fptype var_E_l;

  BasicSimResults() : success( false ) { }

};

#endif // SIMRESULTS_H_INCLUDED
