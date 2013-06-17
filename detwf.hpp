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

#ifndef DETERMINANTAL_WAVEFUNCTIONS_H_INCLUDED
#define DETERMINANTAL_WAVEFUNCTIONS_H_INCLUDED

#include <vector>
#include <memory>

#define EIGEN_NO_AUTOMATIC_RESIZING
#include <eigen3/Eigen/Core>

#include "options.hpp"
#include "macros.h"
#include "fptype.hpp"
#include "lattice.hpp"


class SingleParticleHamiltonian final {

  public:

    // the number of lattice sites
    const unsigned int L;

  private:

    // the single particle Hamiltonian under p.h. transformation
    Eigen::MatrixXfp int_H;

    // masks the vpar operators
    std::vector<Eigen::MatrixXfp> int_V;

  public:

    SingleParticleHamiltonian( unsigned int L_init );

    void add_anyterm(  const Eigen::MatrixXfp& term );
    void add_vparterm( const Eigen::MatrixXfp& mask, fptype vpar );

    const Eigen::MatrixXfp& H() const;
    const std::vector<Eigen::MatrixXfp>& V() const;
};


class DeterminantalWavefunction final {

  private:

    // the single particle Hamiltonian
    const SingleParticleHamiltonian int_spHam;

  public:

    // the number of particles (= number of states in the slater determinant)
    const unsigned int Np;

  private:

    // the eigenstates
    Eigen::MatrixXfp int_U;

    // the associated energies
    Eigen::VectorXfp int_epsilon;

    // resulting matrices for the vpar operators
    std::vector<Eigen::MatrixXfp> int_A;

  public:

    DeterminantalWavefunction(
      const SingleParticleHamiltonian& spHam_init, unsigned int Np_init );

    bool is_openshell() const;

    const SingleParticleHamiltonian& spHam() const;

    const Eigen::MatrixXfp& U() const;
    Eigen::Block<const Eigen::MatrixXfp> M() const;
    const Eigen::VectorXfp& epsilon() const;
    const std::vector<Eigen::MatrixXfp>& A() const;
};


DeterminantalWavefunction build_detwf(
  const std::shared_ptr<Lattice>& lat, unsigned int Ne,
  const std::vector<double>& t,
  const std::vector<double>& Delta, optpairsym_t pairsym,
  double mu, double mu_m
);


double calc_tbdetwf_chempot(
  const std::shared_ptr<Lattice>& lat, unsigned int Ne,
  const std::vector<double>& t
);

#endif // DETERMINANTAL_WAVEFUNCTIONS_H_INCLUDED
