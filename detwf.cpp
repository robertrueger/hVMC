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


VariationalHamiltonian::VariationalHamiltonian( unsigned int L_init )
  : L( L_init ), int_H( Eigen::MatrixXfp::Zero( 2 * L, 2 * L ) ) { }


void VariationalHamiltonian::add_anyterm(  const Eigen::MatrixXfp& term )
{
  int_H += term;
}


void VariationalHamiltonian::add_vparterm(
  const Eigen::MatrixXfp& mask, fptype vpar )
{
  int_H += vpar * mask;
  int_V.push_back( mask );
}


const Eigen::MatrixXfp& VariationalHamiltonian::H() const
{
  return int_H;
}


const vector<Eigen::MatrixXfp>& VariationalHamiltonian::V() const
{
  return int_V;
}


DeterminantalWavefunction::DeterminantalWavefunction(
  const VariationalHamiltonian& varHam_init, unsigned int Np_init )
  : int_varHam( varHam_init ), Np( Np_init ),
    int_U( 2 * int_varHam.L, 2 * int_varHam.L ),
    int_epsilon( 2 * int_varHam.L )
{
  // diagonalize variational Hamiltonian
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      int_varHam.H().cast<double>()
    );
    assert( solver.info() == Eigen::Success );
    int_U = solver.eigenvectors().cast<fptype>();
    int_epsilon = solver.eigenvalues().cast<fptype>();
  }

  // define perturbation theory mask
  Eigen::ArrayXfp ptmask
    = Eigen::ArrayXfp::Zero( 2 * int_varHam.L,  2 * int_varHam.L );
  for ( Lattice::spindex eta = 0; eta < 2 * int_varHam.L; ++eta ) {
    for ( Lattice::spindex nu = 0; nu < 2 * int_varHam.L; ++nu ) {
      if ( eta >= Np && nu < Np ) {
        ptmask( eta, nu ) = 1.f / ( int_epsilon( nu ) - int_epsilon( eta ) );
      } else {
        ptmask( eta, nu ) = 0.f;
      }
    }
  }

  // calculate the A matrices of the variational parameters
  for ( auto it = int_varHam.V().begin(); it != int_varHam.V().end(); ++it ) {
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


const VariationalHamiltonian& DeterminantalWavefunction::varHam() const
{
  return int_varHam;
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
  const std::vector<double>& t,
  const std::vector<double>& Delta, optpairsym_t pairsym,
  double mu, double mu_m )
{
  // check if we have the correct number of variational parameters
  assert( t.size() == 3 );
  assert( Delta.size() == 4 );


  VariationalHamiltonian varHam( lat->L );
  Eigen::MatrixXfp varHam_mask = Eigen::MatrixXfp::Zero( 2 * lat->L, 2 * lat->L );

  // nearest neighbour hopping
  // (NOT a variational parameters, as it determines the energy scale!)
  for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
    vector<Lattice::spindex> l_Xnn = lat->get_Xnn( l, 1 );
    for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
      varHam_mask( l, *it ) =
        ( lat->get_spindex_type( l ) == Lattice::spindex_type::up ) ? -1.f : 1.f;
    }
  }
  varHam.add_anyterm( t[0] * varHam_mask );
  varHam_mask.setZero();

  // 2nd and 3rd nearest neighbor hopping
  for ( unsigned int X = 2; X <= 3; ++X ) {
    for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
      vector<Lattice::spindex> l_Xnn = lat->get_Xnn( l, X );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        varHam_mask( l, *it ) =
          ( lat->get_spindex_type( l ) == Lattice::spindex_type::up ) ? -1.f : 1.f;
      }
    }
    varHam.add_vparterm( varHam_mask, t[X - 1] );
    varHam_mask.setZero();
  }

  // onsite BCS pairing
  for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
    varHam_mask( l, lat->get_linked_spindex( l ) ) = +1.f;
  }
  varHam.add_vparterm( varHam_mask, Delta[0] );
  varHam_mask.setZero();

  // 1st, 2nd and 3rd nearest neighbor BCS pairing
  for ( unsigned int X = 1; X <= 3; ++X ) {
    for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
      vector<Lattice::spindex> l_Xnn = lat->get_Xnn( l, X );
      for ( auto it = l_Xnn.begin(); it != l_Xnn.end(); ++it ) {
        varHam_mask( l, lat->get_linked_spindex( *it ) )
          = +1.f * lat->pairsym_modifier( pairsym, l, *it );
      }
    }
    varHam.add_vparterm( varHam_mask, Delta[X] );
    varHam_mask.setZero();
  }

  // chemical potential
  varHam_mask.diagonal().head( lat->L ).array() -= 1.f;
  varHam_mask.diagonal().tail( lat->L ).array() += 1.f;
  varHam.add_vparterm( varHam_mask, mu );
  varHam_mask.setZero();

  // site and spin dependent chemical potential to introduce magnetism
  for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
    varHam_mask( l, l ) =
      lat->get_index_sublattice( lat->get_index_from_spindex( l ) ) == 0 ?
      -0.5 :
      0.5;
  }
  varHam.add_vparterm( varHam_mask, mu_m );
  varHam_mask.setZero();

  // determine how many particles we have after the p.-h. transformation
  assert( Ne % 2 == 0 );
  const unsigned int Np = Ne / 2 + ( lat->L - Ne / 2 );

  return DeterminantalWavefunction( varHam, Np );
}



double calc_tbdetwf_chempot(
  const std::shared_ptr<Lattice>& lat, unsigned int Ne,
  const std::vector<double>& t )
{
  Eigen::MatrixXfp H_tb_nopht = Eigen::MatrixXfp::Zero( 2 * lat->L, 2 * lat->L );

  for ( Lattice::spindex l = 0; l < 2 * lat->L; ++l ) {
    for ( unsigned int X = 1; X <= 3; ++X ) {
      vector<Lattice::spindex> l_Xnn = lat->get_Xnn( l, X );
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
