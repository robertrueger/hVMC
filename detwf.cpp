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

#include "detwf.hpp"

#include <iostream>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Eigenvalues>

using namespace std;


SingleParticleOrbitals wf_tight_binding(
  const vector<double>& t,
  unsigned int N, const shared_ptr<Lattice>& lat )
{
  // make sure we pass the right number of parameters
  assert( t.size() == 3  );

  Eigen::MatrixXfp H_tb = Eigen::MatrixXfp::Zero( 2*lat->L, 2*lat->L );

  // t hopping
  vector<unsigned int> l_Xnn;
  for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
    for ( unsigned int X = 1; X <= 3; ++X ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        H_tb( l, *it ) += ( l < lat->L ? -1.f : 1.f ) * t[X - 1];
      }
    }
  }

  // chemical potential
  const fptype mu = -0.65; // TODO: calculate (Robert Rueger, 2013-04-19 10:49)
  H_tb.diagonal().head( lat->L ).array() -= mu;
  H_tb.diagonal().tail( lat->L ).array() += mu;

#if VERBOSE >= 1
  cout << "wf_tight_binding() : Hamiltonian in single particle basis ="
       << endl << H_tb << endl;
#endif

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXfp> eigensolver( H_tb );
  assert( eigensolver.info() == Eigen::Success );

  assert( N % 2 == 0 );

  const Eigen::MatrixXfp& M
    = eigensolver.eigenvectors().topLeftCorner( 2 * lat->L, N );

#if VERBOSE >= 1
  cout << "wf_tight_binding() : eigenenergies = " << endl
       << eigensolver.eigenvalues().transpose() << endl
       << "wf_tight_binding() : eigenvectors = " << endl
       << eigensolver.eigenvectors() << endl

       << "wf_tight_binding() : M = " << endl << M << endl
       << "wf_tight_binding() : slater determinant ground state energy = "
       << eigensolver.eigenvalues().head( N ).sum()  << endl;
#endif

  return SingleParticleOrbitals( M, eigensolver.eigenvalues(), false );
}
