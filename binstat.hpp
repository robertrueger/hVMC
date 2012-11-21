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

#ifndef BINSTAT_H_INCLUDED
#define BINSTAT_H_INCLUDED

#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <stdexcept>

#if VERBOSE >= 1
# include <iostream>
#endif

#include "macros.h"
#include "fptype.hpp"


template <typename T>
class Bin final
{

  private:

    std::vector<T> v;

  public:

    Bin( size_t num_binmcs ) : v( std::vector<T>( num_binmcs ) ) { }
    Bin( size_t num_binmcs, T init ) : v( std::vector<T>( num_binmcs, init ) ) { }

    Bin<T>& operator=( const Bin<T>& other ) {
      if ( v.size() == other.v.size() ) {
        v = other.v;
      } else {
        throw std::logic_error("assignment of bins with different sizes");
      }
      return *this;
    }

    size_t size() const {
      return v.size();
    }


    T operator[]( size_t i ) const {
      return v.at( i );
    }
    T& operator[]( size_t i ) {
      return v.at( i );
    }


    typename std::vector<T>::const_iterator begin() const {
      return v.begin();
    }
    typename std::vector<T>::iterator begin() {
      return v.begin();
    }

    typename std::vector<T>::const_iterator end() const {
      return v.end();
    }
    typename std::vector<T>::iterator end() {
      return v.end();
    }


    typename std::vector<T>::const_reverse_iterator rbegin() const {
      return v.rbegin();
    }
    typename std::vector<T>::reverse_iterator rbegin() {
      return v.rbegin();
    }

    typename std::vector<T>::const_reverse_iterator rend() const {
      return v.rend();
    }
    typename std::vector<T>::reverse_iterator rend() {
      return v.rend();
    }
};


template <typename T>
class BinnedData final
{

  private:

    std::vector< Bin<T> > data;

  public:

    BinnedData( size_t num_bins, size_t num_binmcs )
      : data( std::vector< Bin<T> >( num_bins, Bin<T>( num_binmcs ) ) ) {
      if ( num_bins == 0 || num_binmcs == 0 ) {
        throw std::logic_error( "num_bins and num_binmcs must be >0" );
      }
    };
    BinnedData( size_t num_bins, size_t num_binmcs, T init )
      : data( std::vector< Bin<T> >( num_bins, Bin<T>( num_binmcs, init ) ) ) {
      if ( num_bins == 0 || num_binmcs == 0 ) {
        throw std::logic_error( "num_bins and num_binmcs must be >0" );
      }
    };


    size_t num_bins() const {
      return data.size();
    }
    size_t num_binmcs() const {
      return data[0].size();
    }


    const Bin<T>& operator[]( size_t i ) const {
      return data.at( i );
    }
    Bin<T>& operator[]( size_t i ) {
      return data.at( i );
    }


    typename std::vector< Bin<T> >::const_iterator begin() const {
      return data.begin();
    }
    typename std::vector< Bin<T> >::iterator begin() {
      return data.begin();
    }

    typename std::vector< Bin<T> >::const_iterator end() const {
      return data.end();
    }
    typename std::vector< Bin<T> >::iterator end() {
      return data.end();
    }


    typename std::vector< Bin<T> >::const_reverse_iterator rbegin() const {
      return data.rbegin();
    }
    typename std::vector< Bin<T> >::reverse_iterator rbegin() {
      return data.rbegin();
    }

    typename std::vector< Bin<T> >::const_reverse_iterator rend() const {
      return data.rend();
    }
    typename std::vector< Bin<T> >::reverse_iterator rend() {
      return data.rend();
    }
};


struct BinnedDataStatistics final {

  fptype mean;
  fptype sigma_mean;

  fptype variance;

  // TODO: autocorrelation (Robert Rueger, 2012-11-12 14:14)
  // fptype t_autocorr;
  // std::vector<fptype> f_autocorr;

};


template <typename T>
BinnedDataStatistics run_bindat_statanalysis( const BinnedData<T>& bd )
{
  BinnedDataStatistics stat;

  // calculate the mean value of the individual bins
  std::vector<fptype> binmean( bd.num_bins() );

  // TODO: avoid for loop in favor of STL (Robert Rueger, 2012-11-13 11:55)
  for ( unsigned int bin = 0; bin < binmean.size(); ++bin ) {
    assert( bd[bin].size() == bd[0].size() );

#if VERBOSE >= 2
    std::cout << "run_bindat_statanalysis() : bin[" << bin << "] data = "
              << std::endl;
    std::copy( bd[bin].begin(), bd[bin].end(),
               std::ostream_iterator<T>( std::cout, " " ) );
    std::cout << std::endl;
#endif
    
    binmean[bin] =
      static_cast<fptype>( accumulate( bd[bin].begin(), bd[bin].end(), 0.f ) ) /
      static_cast<fptype>( bd[bin].size() );
  }

#if VERBOSE >= 2
  std::cout << "run_bindat_statanalysis() : binmeans are " << std::endl;
  std::copy ( binmean.begin(), binmean.end(),
              std::ostream_iterator<fptype>( std::cout, " " ) );
  std::cout << std::endl;
#endif

  // calculate the average of the bins' mean values
  stat.mean = accumulate( binmean.begin(), binmean.end(), 0.f ) /
              static_cast<fptype>( binmean.size() );

  // calculate the variance of the observable
  stat.variance =
    static_cast<fptype>( binmean.size() )
    / static_cast<fptype>( binmean.size() - 1 )
    * (
      accumulate( binmean.begin(), binmean.end(), 0.f,
  []( fptype sum, fptype m ) {
    return sum + m * m;
  } )
      / static_cast<fptype>( binmean.size() )
      - stat.mean * stat.mean
    );

  // uncertainty of the mean is sqrt(variance / bins)
  stat.sigma_mean = sqrt( stat.variance / static_cast<fptype>( binmean.size() ) );

  // TODO: autocorrelation (Robert Rueger, 2012-11-01 16:28)

  return stat;
}

#endif // BINSTAT_H_INCLUDED
