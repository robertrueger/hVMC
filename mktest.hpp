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

#ifndef MANN_KENDALL_TEST_H_INCLUDED
#define MANN_KENDALL_TEST_H_INCLUDED

#include <deque>
#include <cmath>

#include <boost/math/special_functions/sign.hpp>


template <typename T>
class MannKendall final {

  private:

    std::deque<T> data;
    long S;

  public:

    MannKendall() : S( 0 ) { }

    size_t size() const {
      return data.size();
    }

    void push_back( T x ) {
      for ( auto it = data.begin(); it != data.end(); ++it ) {
        S += boost::math::sign( x - *it );
      }
      data.push_back( x );
    }

    void remove_front() {
      T deleted = data.front();
      data.pop_front();
      for ( auto it = data.begin(); it != data.end(); ++it ) {
        S -= boost::math::sign( *it - deleted );
      }
    }

    void clear() {
      data.clear();
      S = 0;
    }

    double test() {
      return std::abs(
        static_cast<double>( S ) /
        static_cast<double>( data.size() * ( data.size() - 1 ) / 2 )
      );
    }
};

#endif // MANN_KENDALL_TEST_H_INCLUDED
