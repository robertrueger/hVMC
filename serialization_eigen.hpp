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

#ifndef SERIALIZATION_EIGEN_H_INCLUDED
#define SERIALIZATION_EIGEN_H_INCLUDED

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

// adapted from the excellent stackoverflow.com answer by Jakob
// Profile: http://stackoverflow.com/users/672634/jakob
//  Answer: http://stackoverflow.com/questions/12580579#12618789

namespace boost
{
  namespace serialization
  {

    template<
      class Archive,
      class _Scalar, int _Rows, int _Cols,
      int _Options, int _MaxRows, int _MaxCols
    >
    void serialize(
      Archive& ar,
      Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& t,
      const unsigned int
    )
    {
      int rows = t.rows(), cols = t.cols();
      ar & rows;
      ar & cols;
      if ( rows * cols != t.size() ) {
        t.resize( rows, cols );
      }
      ar & make_array(t.data(), t.size());
    }

  }
}

#endif // SERIALIZATION_EIGEN_H_INCLUDED
