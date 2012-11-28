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

#include "detwf.hpp"
using namespace std;


Eigen::MatrixXfp wf_nntb(
    const vector<fptype>& t,
    unsigned int N, Lattice* const lat )
{
  Eigen::MatrixXfp H_tb_nospin = Eigen::MatrixXfp::Zero( lat->L, lat->L );

  vector<unsigned int> l_Xnn;
  for ( unsigned int l = 0; l < lat->L; ++l ) {
    for ( unsigned int X = 1; X <= t.size(); ++X ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        H_tb_nospin( l, *it ) += -1.f * t[X - 1];
      }
    }
  }

#if VERBOSE >= 1
  cout << "wf_nntb() : spinless Hamiltonian in single particle basis =" << endl;
  cout << H_tb_nospin << endl;
#endif

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXfp> eigensolver( H_tb_nospin );
  assert( eigensolver.info() == Eigen::Success );

  Eigen::MatrixXfp M = Eigen::MatrixXfp::Zero( 2 * lat->L, N );

  assert( N % 2 == 0 );
  for ( unsigned int k = 0; k < N / 2; ++k ) {
    M.col( 2 * k ).head( lat->L )     = eigensolver.eigenvectors().col( k );
    M.col( 2 * k + 1 ).tail( lat->L ) = eigensolver.eigenvectors().col( k );
  }

#if VERBOSE >= 1
  cout << "wf_nntb() : M = " << endl << M << endl;
  cout << "wf_nntb() : slater determinant ground state energy = "
       << 2.f * eigensolver.eigenvalues().head( N / 2 ).sum()  << endl;
#endif

  return M;
}

