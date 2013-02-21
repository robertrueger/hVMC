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

#include <iostream>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Eigenvalues>

using namespace std;


SingleParticleOrbitals wf_tight_binding(
  const vector<fptype>& t,
  unsigned int N, const shared_ptr<Lattice>& lat )
{
  // make sure we pass the right number of variational parameters
  assert( t.size() == 3  );

  Eigen::MatrixXfp H_tb_nospin = Eigen::MatrixXfp::Zero( lat->L, lat->L );

  vector<unsigned int> l_Xnn;
  for ( unsigned int l = 0; l < lat->L; ++l ) {
    for ( unsigned int X = 1; X <= 3; ++X ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        H_tb_nospin( l, *it ) += -1.f * t[X - 1];
      }
    }
  }

#if VERBOSE >= 1
  cout << "wf_tight_binding() : spinless Hamiltonian in single particle basis ="
       << endl << H_tb_nospin << endl;
#endif

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXfp> eigensolver( H_tb_nospin );
  assert( eigensolver.info() == Eigen::Success );

  assert( N % 2 == 0 );

/*
  Eigen::MatrixXfp M = Eigen::MatrixXfp::Zero( 2 * lat->L, N );
  for ( unsigned int k = 0; k < N / 2; ++k ) {
    M.col( 2 * k ).head( lat->L )     = eigensolver.eigenvectors().col( k );
    M.col( 2 * k + 1 ).tail( lat->L ) = eigensolver.eigenvectors().col( k );
  }
  Eigen::VectorXfp E = Eigen::VectorXfp::Zero( 2 * lat->L );
  for ( unsigned int k = 0; k < lat->L; ++k ) {
    E( 2 * k )     = eigensolver.eigenvalues()( k );
    E( 2 * k + 1 ) = eigensolver.eigenvalues()( k );
  }
*/

  const Eigen::MatrixXfp& M
    = eigensolver.eigenvectors().topLeftCorner( lat->L, N / 2 );

#if VERBOSE >= 1
  cout << "wf_tight_binding() : M = " << endl << M << endl
       << "wf_tight_binding() : slater determinant ground state energy = "
       << 2.f* eigensolver.eigenvalues().head( N / 2 ).sum()  << endl;
#endif

  //return SingleParticleOrbitals( M, E, false );
  return SingleParticleOrbitals( M, eigensolver.eigenvalues(), true );
}
