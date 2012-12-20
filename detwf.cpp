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


SingleParticleOrbitals wf_tight_binding(
  const vector<cl_fptype>& t,
  cl_uint N, Lattice* lat )
{
  cout << "      > Constructing TB Hamiltonian" << endl;

  Eigen::MatrixXfp H_tb_nospin = Eigen::MatrixXfp::Zero( lat->L, lat->L );

  vector<cl_uint> l_Xnn;
  for ( cl_uint l = 0; l < lat->L; ++l ) {
    for ( cl_uint X = 1; X <= t.size(); ++X ) {
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

  cout << "      > Diagonalizing TB Hamiltonian" << endl;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXfp> eigensolver( H_tb_nospin );
  assert( eigensolver.info() == Eigen::Success );

  assert( N % 2 == 0 );

  /*
    Eigen::MatrixXfp M = Eigen::MatrixXfp::Zero( 2 * lat->L, N );
    for ( cl_uint k = 0; k < N / 2; ++k ) {
      M.col( 2 * k ).head( lat->L )     = eigensolver.eigenvectors().col( k );
      M.col( 2 * k + 1 ).tail( lat->L ) = eigensolver.eigenvectors().col( k );
    }
  */

  const Eigen::MatrixXfp& M
    = eigensolver.eigenvectors().topLeftCorner( lat->L, N / 2 );
#if VERBOSE >= 1
  cout << "wf_tight_binding() : M = " << endl << M << endl
       << "wf_tight_binding() : slater determinant ground state energy = "
       << 2.f* eigensolver.eigenvalues().head( N / 2 ).sum()  << endl;
#endif

  check_openshell( eigensolver.eigenvalues(), N / 2 );   

  return SingleParticleOrbitals( M, true );
}



bool check_openshell( const Eigen::VectorXfp& E, cl_uint N )
{
  cout << "      > Checking for open shell" << endl;

#if VERBOSE >= 1
  cout << "check_openshell() : single particle orbital energies =" << endl;
  for ( cl_uint n = 0; n < E.size(); ++n ) {
    cout << E( n );
    if ( n == N - 1 ) {
      cout << " <-- E_fermi";
    }
    cout << endl;
  }
#endif

  if ( E( N ) - E( N - 1 ) < 0.00001 ) {
    cout << endl;
    cout << "        ERROR: Open shell detected!" << endl;
    cout << "        E_fermi = " << E( N - 1 ) << endl;
    cout << "        Orbital below = " << E( N - 2 ) << endl;
    cout << "        Orbital above = " << E( N ) << endl;
    exit( 1 );
  }

  return true;
}
