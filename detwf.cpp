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

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Eigenvalues>

using namespace std;


SingleParticleHamiltonian::SingleParticleHamiltonian( unsigned int L_init )
  : L( L_init ), int_H( Eigen::MatrixXfp::Zero( 2 * L, 2 * L ) ) { }


void SingleParticleHamiltonian::add_anyterm(  const Eigen::MatrixXfp& term )
{
  int_H += term;
}


void SingleParticleHamiltonian::add_vparterm(
  const Eigen::MatrixXfp& mask, fptype vpar )
{
  int_H += vpar * mask;
  int_V.push_back( mask );
}


const Eigen::MatrixXfp& SingleParticleHamiltonian::H() const
{
  return int_H;
}


const vector<Eigen::MatrixXfp>& SingleParticleHamiltonian::V() const
{
  return int_V;
}


DeterminantalWavefunction::DeterminantalWavefunction(
  const SingleParticleHamiltonian& spHam_init, unsigned int Np_init )
  : int_spHam( spHam_init ), Np( Np_init ),
    int_U( 2 * int_spHam.L, 2 * int_spHam.L ),
    int_epsilon( 2 * int_spHam.L )
{
  // diagonalize single particle Hamiltonian
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      int_spHam.H().cast<double>()
    );
    assert( solver.info() == Eigen::Success );
    int_U = solver.eigenvectors().cast<fptype>();
    int_epsilon = solver.eigenvalues().cast<fptype>();
  }

  // define perturbation theory mask
  Eigen::ArrayXfp ptmask
    = Eigen::ArrayXfp::Zero( 2 * int_spHam.L,  2 * int_spHam.L );
  for ( unsigned int eta = 0; eta < 2 * int_spHam.L; ++eta ) {
    for ( unsigned int nu = 0; nu < 2 * int_spHam.L; ++nu ) {
      if ( eta >= Np && nu < Np ) {
        ptmask( eta, nu ) = 1.f / ( int_epsilon( nu ) - int_epsilon( eta ) );
      } else {
        ptmask( eta, nu ) = 0.f;
      }
    }
  }

  // calculate the A matrices of the variational parameters
  for ( auto it = int_spHam.V().begin(); it != int_spHam.V().end(); ++it ) {
    int_A.push_back(
      int_U *
      ( ( int_U.adjoint() * *it * int_U  ).array() * ptmask ).matrix()
      * int_U.adjoint()
    );
  }
}


bool DeterminantalWavefunction::is_openshell() const
{
  return ( int_epsilon( Np ) - int_epsilon( Np - 1 ) < 0.0001 );
}


const SingleParticleHamiltonian& DeterminantalWavefunction::spHam() const
{
  return int_spHam;
}


const Eigen::MatrixXfp& DeterminantalWavefunction::U() const
{
  return int_U;
}


Eigen::Block<const Eigen::MatrixXfp> DeterminantalWavefunction::M() const
{
  // TODO: return ColsBlockXpr ??? (Robert Rueger, 2013-04-25 14:20)
  return int_U.topLeftCorner( int_U.rows(), Np );
}


const Eigen::VectorXfp& DeterminantalWavefunction::epsilon() const
{
  return int_epsilon;
}


const vector<Eigen::MatrixXfp>& DeterminantalWavefunction::A() const
{
  return int_A;
}



DeterminantalWavefunction build_detwf(
  const std::shared_ptr<Lattice>& lat, unsigned int Ne,
  const std::vector<double>& t, const std::vector<double>& Delta, double mu )
{
  // check if we have the correct number of variational parameters
  assert( t.size() == 3 );
  assert( Delta.size() == 4 );


  SingleParticleHamiltonian spHam( lat->L );
  Eigen::MatrixXfp spHam_mask = Eigen::MatrixXfp::Zero( 2 * lat->L, 2 * lat->L );
  vector<unsigned int> l_Xnn;

  // nearest neighbour hopping
  // (NOT a variational parameters, as it determines the energy scale!)
  for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
    lat->get_Xnn( l, 1, &l_Xnn );
    for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
      spHam_mask( l, *it ) = ( l < lat->L ? -1.f : 1.f );
    }
  }
  spHam.add_anyterm( t[0] * spHam_mask );
  spHam_mask.setZero();

  // 2nd and 3rd nearest neighbor hopping
  for ( unsigned int X = 2; X <= 3; ++X ) {
    for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        spHam_mask( l, *it ) = ( l < lat->L ? -1.f : 1.f );
      }
    }
    spHam.add_vparterm( spHam_mask, t[X - 1] );
    spHam_mask.setZero();
  }

  // onsite BCS pairing
  for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
    if ( l < lat->L ) {
      // l is a spin up site
      spHam_mask( l, l + lat->L ) = +1.f;
    } else {
      // l is a spin down site
      spHam_mask( l, l - lat->L ) = +1.f;
    }
  }
  spHam.add_vparterm( spHam_mask, Delta[0] );
  spHam_mask.setZero();

  // 1st, 2nd and 3rd nearest neighbor BCS pairing
  for ( unsigned int X = 1; X <= 3; ++X ) {
    for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        if ( l < lat->L ) {
          // l is a spin up site
          spHam_mask( l, *it + lat->L ) = +1.f;
        } else {
          // l is a spin down site
          spHam_mask( l, *it - lat->L ) = +1.f;
        }
      }
    }
    spHam.add_vparterm( spHam_mask, Delta[X] );
    spHam_mask.setZero();
  }

  // chemical potential
  spHam_mask.diagonal().head( lat->L ).array() -= 1.f;
  spHam_mask.diagonal().tail( lat->L ).array() += 1.f;
  spHam.add_vparterm( spHam_mask, mu );
  spHam_mask.setZero();


  // determine how many particles we have after the p.-h. transformation
  assert( Ne % 2 == 0 );
  const unsigned int Np = Ne / 2 + ( lat->L - Ne / 2 );

  return DeterminantalWavefunction( spHam, Np );
}



double calc_tbdetwf_chempot(
  const std::shared_ptr<Lattice>& lat, unsigned int Ne,
  const std::vector<double>& t )
{
  Eigen::MatrixXfp H_tb_nopht = Eigen::MatrixXfp::Zero( 2 * lat->L, 2 * lat->L );

  vector<unsigned int> l_Xnn;
  for ( unsigned int l = 0; l < 2 * lat->L; ++l ) {
    for ( unsigned int X = 1; X <= 3; ++X ) {
      lat->get_Xnn( l, X, &l_Xnn );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        H_tb_nopht( l, *it ) -= t[X - 1];
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXfp> H_tb_nopht_solver( H_tb_nopht );
  assert( H_tb_nopht_solver.info() == Eigen::Success );

  return 0.5 * ( H_tb_nopht_solver.eigenvalues()( Ne )
                 + H_tb_nopht_solver.eigenvalues()( Ne - 1 ) );
}
