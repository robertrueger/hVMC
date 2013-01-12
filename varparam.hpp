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

#ifndef VARPARAM_H_INCLUDED
#define VARPARAM_H_INCLUDED

#include <vector>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include "macros.h"
#include "fptype.hpp"
#include "options.hpp"


struct VariationalParameters final {

    std::vector<fptype> determinantal;
    std::vector<fptype> jastrow;

  private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive& ar, const unsigned int ) {
      ar & determinantal;
      ar & jastrow;
    }
};

VariationalParameters get_initial_varparam( const Options& opts );

#endif // VARPARAM_H_INCLUDED
