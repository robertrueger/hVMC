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

#include "obs_dkdkp.hpp"

#if VERBOSE >= 1
# include <iostream>
#endif

#include <boost/mpi/collectives.hpp>

#include "macros.h"
#include "serialization_eigen.hpp"

using namespace std;
namespace mpi = boost::mpi;


ObservableDeltaKDeltaKPrime::ObservableDeltaKDeltaKPrime()
  : Observable( OBSERVABLE_DELTAK_DELTAKPRIME ),
    this_bin_num_measurements( 0 ) { }


void ObservableDeltaKDeltaKPrime::measure( HubbardModelVMC& model )
{
  const Eigen::VectorXfp& Dk_current = model.Delta_k();

  if ( DkDkp_sum.size() == 0 ) {
    // first use of DkDkp_sum
    DkDkp_sum.setZero( Dk_current.size(), Dk_current.size() );
  } else {
    assert( DkDkp_sum.cols() == Dk_current.size() );
    assert( DkDkp_sum.rows() == Dk_current.size() );
  }

  DkDkp_sum += Dk_current * Dk_current.transpose();
  ++this_bin_num_measurements;

#if VERBOSE >= 1
  cout << "ObservableDeltaKDeltaKPrime::measure() : DkDkp_sum = " << endl
       << DkDkp_sum << endl;
#endif
}


void ObservableDeltaKDeltaKPrime::completebin()
{
  DkDkp_binmeans.push_back(
    DkDkp_sum / static_cast<fptype>( this_bin_num_measurements )
  );

#if VERBOSE >= 1
  cout << "ObservableDeltaKDeltaKPrime::completebin() : binmean = " << endl
       << DkDkp_binmeans.rbegin()->transpose() << endl;
#endif

  DkDkp_sum.setZero();
  this_bin_num_measurements = 0;
}


void ObservableDeltaKDeltaKPrime::collect_and_write_results(
  const mpi::communicator& mpicomm,
  MCCResults& results ) const
{
  assert( mpicomm.rank() == 0 );
  vector< vector<Eigen::MatrixXfp> > binmeans_collector;
  mpi::gather( mpicomm, DkDkp_binmeans, binmeans_collector, 0 );

  vector< Eigen::MatrixXfp > DkDkp_binmeans_all;
  for ( auto it = binmeans_collector.begin();
        it != binmeans_collector.end();
        ++it ) {
    DkDkp_binmeans_all.insert( DkDkp_binmeans_all.end(), it->begin(), it->end() );
  }
  assert( !DkDkp_binmeans_all.empty() );


  Eigen::MatrixXfp DkDkp_binmeans_all_sum;
  DkDkp_binmeans_all_sum.setZero(
    DkDkp_binmeans_all[0].rows(),
    DkDkp_binmeans_all[0].cols()
  );
  for ( auto it = DkDkp_binmeans_all.begin();
        it != DkDkp_binmeans_all.end();
        ++it ) {
    DkDkp_binmeans_all_sum += *it;
#if VERBOSE >= 2
    cout << it->transpose() << endl;
#endif
  }
  results.Deltak_Deltakprime
    = DkDkp_binmeans_all_sum / static_cast<fptype>( DkDkp_binmeans_all.size() );
}


void ObservableDeltaKDeltaKPrime::send_results_to_master(
  const mpi::communicator& mpicomm ) const
{
  assert( mpicomm.rank() != 0 );
  mpi::gather( mpicomm, DkDkp_binmeans, 0 );
}
