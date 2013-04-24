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
  const vector<double>& t_ext,
  unsigned int N, const shared_ptr<Lattice>& lat,
  bool is_master )
{
  // TODO: remove and use t_ext (Robert Rueger, 2013-04-22 20:28)
  vector<double> t = { 1.0, -0.3862688, 0 };

  // make sure we pass the right number of parameters
  assert( t.size() == 3  );

  // diagonalize single particle Hamiltonian under p.-h. transformation

  Eigen::MatrixXd H_tb = Eigen::MatrixXd::Zero( 2 * lat->L, 2 * lat->L );

  // t hopping
  vector<unsigned int> l_Xnn;
  for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
    for ( unsigned int X = 1; X <= 3; ++X ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        H_tb( l, *it ) += ( l < lat->L ? -1.0 : 1.0 ) * t[X - 1];
      }
    }
  }

  // BCS term
  const vector<double> Delta = { -0.1800535, 0.3728459, 0.0, 0.0 };
  // TODO: make an option (Robert Rueger, 2013-04-22 16:18)
  // TODO: make optimizable (Robert Rueger, 2013-04-22 16:18)
  for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {

    // onsite BCS term
    if ( l < lat->L ) {
      // l is a spin up site
      H_tb( l, l + lat->L ) += Delta[ 0 ];
    } else {
      // l is a spin down site
      H_tb( l, l - lat->L ) += Delta[ 0 ];
    }

    // neighboring sites BCS term
    for ( unsigned int X = 1; X <= 3; ++X ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        if ( l < lat->L ) {
          // l is a spin up site
          H_tb( l, *it + lat->L ) += Delta[ X ];
        } else {
          // l is a spin down site
          H_tb( l, *it - lat->L ) += Delta[ X ];
        }
      }
    }
  }

  // chemical potential
  const double mu = 2.017839;
  // TODO: make optimizable (Robert Rueger, 2013-04-22 16:18)
  H_tb.diagonal().head( lat->L ).array() -= mu;
  H_tb.diagonal().tail( lat->L ).array() += mu;

#if VERBOSE >= 1
  cout << "wf_tight_binding() : Hamiltonian in single particle basis ="
       << endl << H_tb << endl;
#endif

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> H_tb_solver( H_tb );
  assert( H_tb_solver.info() == Eigen::Success );

  // determine how many particles we have after the p.-h. transformation
  assert( N % 2 == 0 );
  const unsigned int Np = N / 2 + ( lat->L - N / 2 );

  const Eigen::MatrixXd& M
    = H_tb_solver.eigenvectors().topLeftCorner( 2 * lat->L, Np );

#if VERBOSE >= 1
  cout << "wf_tight_binding() : eigenenergies = " << endl
       << H_tb_solver.eigenvalues().transpose() << endl
       << "wf_tight_binding() : eigenvectors = " << endl
       << H_tb_solver.eigenvectors() << endl

       << "wf_tight_binding() : M = " << endl << M << endl
       << "wf_tight_binding() : slater determinant ground state energy = "
       << H_tb_solver.eigenvalues().head( N ).sum()  << endl;
#endif

  // ----- 3.: check for open shell
  if ( H_tb_solver.eigenvalues()( Np ) - H_tb_solver.eigenvalues()( Np - 1 ) < 0.00001 ) {
    if ( is_master ) {
      cout << endl;
      cout << "   ERROR: Open shell detected!" << endl;
      cout << "     E_fermi = " << H_tb_solver.eigenvalues()( Np - 1 ) << endl;
      cout << "     Orbital below = " << H_tb_solver.eigenvalues()( Np - 2 ) << endl;
      cout << "     Orbital above = " << H_tb_solver.eigenvalues()( Np ) << endl;
      cout << endl;
    }
    exit( 1 );
  }

  return SingleParticleOrbitals(
           M.cast<fptype>(),
           H_tb_solver.eigenvalues().cast<fptype>()
         );
}
