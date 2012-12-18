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

#include "hmodvmc_cl.hpp"
using namespace std;



HubbardModelVMC_CL::HubbardModelVMC_CL(
  const cl::Context& clCtx_init,
  mt19937 rng_init,
  Lattice* const lat_init,
  const SingleParticleOrbitals& M_init,
  const Jastrow& v_init,
  cl_uint N_init,
  cl_uint update_hop_maxdist_init,
  const vector<cl_fptype>& t_init,
  cl_fptype U_init,
  cl_fptype W_deviation_target_init,
  cl_uint updates_until_W_recalc_init,
  cl_fptype T_deviation_target_init,
  cl_uint updates_until_T_recalc_init )
  : clCtx( clCtx_init ),
    rng( rng_init ),
    lat( lat_init ), M( M_init ), v( v_init ),
    update_hop_maxdist( update_hop_maxdist_init ),
    t( t_init ), U( U_init ),
    econf( ElectronConfiguration( lat, N_init, &rng ) ),
    Wbu(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, N_init )
    ),
    Wd(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
   Wbu_fromdev(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat->L, N_init )
    ),
    Wd_fromdev(
      M_init.ssym ?
      Eigen::MatrixXfp( lat->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),   devWbu_active( nullptr ), devWbu_inactive( nullptr ),
    devWd_active( nullptr ), devWd_inactive( nullptr ),
    T( Eigen::VectorXfp( lat->L ) ),
    completed_mcsteps( 0 ),
    updates_until_W_recalc( updates_until_W_recalc_init ),
    updates_until_T_recalc( updates_until_T_recalc_init ),
    updates_since_W_recalc( 0 ), updates_since_T_recalc( 0 ),
    W_devstat( FPDevStat( W_deviation_target_init ) ),
    T_devstat( FPDevStat( T_deviation_target_init ) )
{
  // setup up everything OpenCL related (queue, program, buffers, etc ...)
  ocl_setup();

  do { // loop to get a large overlap

    bool invertible;

    do { // loop to get any overlap at all

      econf.distribute_random();

      if ( M.ssym == true ) {

        Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Du().transpose() );
        invertible = ( lu_decomp.rank() == lu_decomp.rows() );
        if ( invertible ) {
          Wbu.noalias() = lu_decomp.solve( M.orbitals.transpose() ).transpose();
        } else {
          continue;
        }

        lu_decomp.compute( calc_Dd().transpose() );
        invertible &= ( lu_decomp.rank() == lu_decomp.rows() );
        if ( invertible ) {
          Wd.noalias() = lu_decomp.solve( M.orbitals.transpose() ).transpose();
        }

      } else {

        Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Db() );
        invertible = ( lu_decomp.rank() == lu_decomp.rows() );
        if ( invertible ) {
          Wbu.noalias() = lu_decomp.solve( M.orbitals.transpose() ).transpose();
        }

      }

#if VERBOSE >= 1
      cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : "
           << "matrix D is not invertible!" << endl;
#endif

    } while ( !invertible );
    // initialize the electrons so that D is invertible
    // (there must be a non-zero overlap between the slater det and |x>)

    // calculate the vector T from scratch
    T = calc_new_T();

    // repeat everything if the initial state has a very low overlap
    // (it's bad for floating point precision before the first recalc)
  } while (

    Wbu.array().square().sum()
    / static_cast<cl_fptype>( Wbu.size() ) > 10.f ||

    (
      M.ssym == true ?

      Wd.array().square().sum()
      / static_cast<cl_fptype>( Wd.size() ) > 10.f :

      false
    ) ||

    T.array().square().sum() / static_cast<cl_fptype>( T.size() ) > 10.f

  );

#if VERBOSE >= 1
  cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : "
       << "calculated initial matrix W =" << endl << Wbu << endl;
  if ( M.ssym == true ) {
    cout << "----->" << endl << Wd << endl;
  }
#endif

  // enqueue transfer of Wbu and Wd to the device
  clQ_Wbu.enqueueWriteBuffer(
    *devWbu_active, CL_FALSE,
    0, Wbu.size() * sizeof( cl_fptype ), Wbu.data()
  );
  if ( M.ssym == true ) {
    clQ_Wd.enqueueWriteBuffer(
      *devWd_active, CL_FALSE,
      0, Wd.size() * sizeof( cl_fptype ), Wd.data()
    );
  }
}



HubbardModelVMC_CL::~HubbardModelVMC_CL()
{
  delete lat;
}



void HubbardModelVMC_CL::ocl_setup()
{
  // read the programs source code from hmodvmc_cl.cl
  string hmodvmc_sourcecode =
    get_file_contents( get_hVMC_dir() / "hmodvmc_cl.cl" );
  cl::Program::Sources hmodvmc_clsources;
  hmodvmc_clsources.push_back(
    make_pair( hmodvmc_sourcecode.c_str(), hmodvmc_sourcecode.length() + 1 )
  );

  // construct the program
  clPrg = cl::Program( clCtx, hmodvmc_clsources );

  // extract the devices from the context and make sure there is only one
  vector<cl::Device> used_devices = clCtx.getInfo<CL_CONTEXT_DEVICES>();
  assert( used_devices.size() == 1 );

  // build the program
  try {
    string compileropts = "";
#ifdef USE_FP_DBLPREC
    compileropts += "-DUSE_FP_DBLPREC_OPENCL";
#endif
    clPrg.build( used_devices, compileropts.c_str() );
  } catch ( cl::Error& e ) {
    cerr << "Build failed! " << e.what() << " (" << e.err() << ")" << endl;
    cerr << "Build log:" << endl;
    cerr << clPrg.getBuildInfo<CL_PROGRAM_BUILD_LOG>( used_devices[0] ) << endl;
    throw e;
  }

  // setup the command queue
  clQ_Wbu = cl::CommandQueue(
              clCtx, used_devices[0]
            );
  clQ_Wd = cl::CommandQueue(
             clCtx, used_devices[0]
           );

  // set up W buffers on the device
  devWbu_1 = cl::Buffer(
               clCtx, CL_MEM_READ_WRITE,
               Wbu.size() * sizeof( cl_fptype )
             );
  devWbu_2 = cl::Buffer(
               clCtx, CL_MEM_READ_WRITE,
               Wbu.size() * sizeof( cl_fptype )
             );
  devWbu_active   = &devWbu_1;
  devWbu_inactive = &devWbu_2;

  if ( M.ssym == true ) {
    devWd_1 = cl::Buffer(
                clCtx, CL_MEM_READ_WRITE,
                Wd.size() * sizeof( cl_fptype )
              );
    devWd_2 = cl::Buffer(
                clCtx, CL_MEM_READ_WRITE,
                Wd.size() * sizeof( cl_fptype )
              );
    devWd_active   = &devWd_1;
    devWd_inactive = &devWd_2;
  }

  // set up Kernel objects
  clK_update_devW = cl::Kernel( clPrg, "update_W" );
}



void HubbardModelVMC_CL::mcs()
{
#if VERBOSE >= 1
  cout << "HubbardModelVMC_CL::mcs() : starting new Monte Carlo step!" << endl;
#endif

  // perform a number of metropolis steps equal to the number of electrons
  for ( cl_uint s = 0; s < econf.N(); ++s ) {
#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::mcs() : "
         << "Monte Carlo step = " << completed_mcsteps
         << ", Metropolis step = " << s << endl;
#endif
    metstep();
  }
  ++completed_mcsteps;
}



void HubbardModelVMC_CL::equilibrate( cl_uint N_mcs_equil )
{
  for ( cl_uint n = 0; n < N_mcs_equil; ++n ) {
    mcs();
  }
  completed_mcsteps -= N_mcs_equil;
}



bool HubbardModelVMC_CL::metstep()
{
  // let the electron configuration propose a random hop
  const ElectronHop& phop = econf.propose_random_hop( update_hop_maxdist );


  // check if the hop is possible (hopto site must be empty)
  if ( phop.possible == false ) {

    // hop is not possible, rejected!
#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::metstep() : hop impossible!" << endl;
#endif
    return false;

  } else { // hop possible!

    // enqueue transfer of W_lk from the device to the host
    cl_fptype W_lk;
    if ( M.ssym == true && phop.k >= econf.N() / 2 ) {
      // the W_lk element is part of devWd_active
      clQ_Wd.enqueueReadBuffer(
        *devWd_active, CL_FALSE,
        ( ( phop.l - lat->L ) * Wd.cols() + ( phop.k - econf.N() / 2 ) )
        * sizeof( cl_fptype ),
        sizeof( cl_fptype ), &W_lk
      );
    } else {
      // the W_lk element is part of devWbu_active
      clQ_Wbu.enqueueReadBuffer(
        *devWbu_active, CL_FALSE,
        ( phop.l * Wbu.cols() + phop.k ) * sizeof( cl_fptype ),
        sizeof( cl_fptype ), &W_lk
      );
    }

    const cl_fptype R_j =
      T( lat->get_spinup_site( phop.l ) )
      / T( lat->get_spinup_site( phop.k_pos ) )
      * v.exp_onsite() / v.exp( phop.l, phop.k_pos );

    // wait for the W_lk transfer to finish
    if ( M.ssym == true && phop.k >= econf.N() / 2 ) {
      clQ_Wd.finish();
    } else {
      clQ_Wbu.finish();
    }

    const cl_fptype accept_prob = R_j * R_j * W_lk * W_lk;

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::metstep() : hop possible -> "
         << "R_j = " << R_j
         << ", sdwf_ratio = " << W_lk
         << ", accept_prob = " << accept_prob << endl;
#endif

    if ( accept_prob >= 1.f ||
         uniform_real_distribution<cl_fptype>( 0.f, 1.f )( rng ) < accept_prob ) {

#if VERBOSE >= 1
      cout << "HubbardModelVMC_CL::metstep() : hop accepted!" << endl;
#endif

      econf.do_hop( phop );

      perform_W_update( phop );
      perform_T_update( phop );

      return true;

    } else { // hop possible but rejected!

#if VERBOSE >= 1
      cout << "HubbardModelVMC_CL::metstep() : hop rejected!" << endl;
#endif

      return false;
    }
  }
}



void HubbardModelVMC_CL::perform_W_update( const ElectronHop& hop )
{
  if ( updates_since_W_recalc >= updates_until_W_recalc ) {

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_W_update() : recalculating W!" << endl;
#endif

    // enqueue update of W on the device
    // (we want to update it to the point of the recalculated W)
    calc_qupdated_W( hop );

    // enqueue transfer of the updated W from the device to the host
    // (waits for the update to finish first)
    clQ_Wbu.enqueueReadBuffer(
      *devWbu_active, CL_FALSE,
      0, Wbu_fromdev.size() * sizeof( cl_fptype ), Wbu_fromdev.data()
    );
    if ( M.ssym == true ) {
      clQ_Wd.enqueueReadBuffer(
        *devWd_active, CL_FALSE,
        0, Wd_fromdev.size() * sizeof( cl_fptype ), Wd_fromdev.data()
      );
    }

    // recalculate W on the host
    calc_new_W();

    // wait for transfers from the device to finish
    clQ_Wbu.finish();
    clQ_Wd.finish();

    // enqueue upload of the recalculated W to the device
    clQ_Wbu.enqueueWriteBuffer(
      *devWbu_active, CL_FALSE,
      0, Wbu.size() * sizeof( cl_fptype ), Wbu.data()
    );
    if ( M.ssym == true ) {
      clQ_Wd.enqueueWriteBuffer(
        *devWd_active, CL_FALSE,
        0, Wd.size() * sizeof( cl_fptype ), Wd.data()
      );
    }

    // compare the recalculated W and the updated W on the host
    cl_fptype dev = calc_deviation( Wbu, Wbu_fromdev );
    if ( M.ssym == true ) {
      dev += calc_deviation( Wd, Wd_fromdev );
    }
    W_devstat.add( dev );

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_W_update() : recalculated W "
         << "with deviation = " << dev << endl;

    if ( dev > W_devstat.target ) {
      cout << "HubbardModelVMC_CL::perform_W_update() : deviation goal for matrix "
           << "W not met!" << endl
           << "HubbardModelVMC_CL::perform_W_update() : approximate W =" << endl
           << Wbu_fromdev << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << Wd_fromdev << endl;
      }
      cout << "HubbardModelVMC_CL::perform_W_update() : exact W =" << endl
           << Wbu << endl;
      if ( M.ssym == true ) {
        cout << "----->" << endl << Wd << endl;
      }
    }
#endif

    assert( dev < W_devstat.target );

    updates_since_W_recalc = 0;

  } else {

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_W_update() : "
         << "performing a quick update of W!" << endl;
#endif

    calc_qupdated_W( hop );
    ++updates_since_W_recalc;
  }
}



Eigen::MatrixXfp HubbardModelVMC_CL::calc_Db() const
{
  assert( M.ssym == false );

  Eigen::MatrixXfp Db( econf.N(), econf.N() );
  for ( cl_uint eid = 0; eid < econf.N(); ++eid ) {
    Db.row( eid ) = M.orbitals.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC_CL::calc_Db() : Db = " << endl << Db << endl;
#endif

  return Db;
}



Eigen::MatrixXfp HubbardModelVMC_CL::calc_Du() const
{
  assert( M.ssym == true );

  Eigen::MatrixXfp Du( econf.N() / 2, econf.N() / 2 );
  for ( cl_uint eid = 0; eid < econf.N() / 2; ++eid ) {
    Du.row( eid ) = M.orbitals.row( econf.get_electron_pos( eid ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC_CL::calc_Du() : Du = " << endl << Du << endl;
#endif

  return Du;
}



Eigen::MatrixXfp HubbardModelVMC_CL::calc_Dd() const
{
  assert( M.ssym == true );

  Eigen::MatrixXfp Dd( econf.N() / 2, econf.N() / 2 );
  for ( cl_uint eid = econf.N() / 2; eid < econf.N(); ++eid ) {
    Dd.row( eid - econf.N() / 2 )
      = M.orbitals.row( lat->get_spinup_site( econf.get_electron_pos( eid ) ) );
  }

#if VERBOSE >= 2
  cout << "HubbardModelVMC_CL::calc_Dd() : Dd = " << endl << Dd << endl;
#endif

  return Dd;
}


void HubbardModelVMC_CL::calc_new_W()
{
  if ( M.ssym == true ) {
    Wbu.noalias() =
      calc_Du().transpose()
      .partialPivLu().solve( M.orbitals.transpose() )
      .transpose();
    Wd.noalias() =
      calc_Dd().transpose()
      .partialPivLu().solve( M.orbitals.transpose() )
      .transpose();
  } else {
    Wbu.noalias() =
      calc_Db().transpose()
      .partialPivLu().solve( M.orbitals.transpose() )
      .transpose();
  }
}



void HubbardModelVMC_CL::calc_qupdated_W( const ElectronHop& hop )
{
  // enqueue the kernel that updates W on the device
  if ( M.ssym == true && hop.k >= econf.N() / 2 ) {
    // devWd must be updated
    clK_update_devW.setArg( 0, *devWd_active );
    clK_update_devW.setArg( 1, *devWd_inactive );
    clK_update_devW.setArg( 2, static_cast<cl_uint>( Wd.cols() ) );
    clK_update_devW.setArg( 3, hop.k - econf.N() / 2 );
    clK_update_devW.setArg( 4, hop.l - lat->L );
    clK_update_devW.setArg( 5, hop.k_pos - lat->L );
    cl::NDRange global( Wd.size() );
    clQ_Wd.enqueueNDRangeKernel(
      clK_update_devW,
      cl::NullRange, global, cl::NullRange
    );
    swap( devWd_active, devWd_inactive );
  } else {
    // devWbu must be updated
    clK_update_devW.setArg( 0, *devWbu_active );
    clK_update_devW.setArg( 1, *devWbu_inactive );
    clK_update_devW.setArg( 2, static_cast<cl_uint>( Wbu.cols() ) );
    clK_update_devW.setArg( 3, hop.k );
    clK_update_devW.setArg( 4, hop.l );
    clK_update_devW.setArg( 5, hop.k_pos );
    cl::NDRange global( Wbu.size() );
    clQ_Wbu.enqueueNDRangeKernel(
      clK_update_devW,
      cl::NullRange, global, cl::NullRange
    );
    swap( devWbu_active, devWbu_inactive );
  }
}



void HubbardModelVMC_CL::perform_T_update( const ElectronHop& hop )
{
  if ( updates_since_T_recalc >= updates_until_T_recalc ) {

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_T_update() : recalculating T!" << endl;
#endif

    updates_since_T_recalc = 0;

    const Eigen::MatrixXfp& T_approx = calc_qupdated_T( hop );
    T = calc_new_T();

    cl_fptype dev = calc_deviation( T_approx, T );
    T_devstat.add( dev );

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_T_update() : recalculated T "
         << "with deviation = " << dev << endl;

    if ( dev > T_devstat.target ) {
      cout << "HubbardModelVMC_CL::perform_T_update() : deviation goal for matrix "
           << "T not met!" << endl
           << "HubbardModelVMC_CL::perform_T_update() : approximate T =" << endl
           << T_approx.transpose() << endl
           << "HubbardModelVMC_CL::perform_T_update() : exact T =" << endl
           << T.transpose() << endl;
    }
#endif

    assert( dev < T_devstat.target );

  } else {

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_T_update() : "
         << "performing a quick update of T!" << endl;
#endif

    ++updates_since_T_recalc;

    T = calc_qupdated_T( hop );

#ifndef NDEBUG
    const Eigen::MatrixXfp& T_chk = calc_new_T();
    cl_fptype dev = calc_deviation( T, T_chk );

# if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_T_update() : "
         << "[DEBUG CHECK] deviation after quick update = " << dev << endl;

    if ( dev > T_devstat.target ) {
      cout << "HubbardModelVMC_CL::perform_T_update() : deviation goal for matrix "
           << "T not met!" << endl
           << "HubbardModelVMC_CL::perform_T_update() : quickly updated T =" << endl
           << T.transpose() << endl
           << "HubbardModelVMC_CL::perform_T_update() : exact T =" << endl
           << T_chk.transpose() << endl;
    }
# endif
#endif

    assert( dev < T_devstat.target );
  }
}



Eigen::VectorXfp HubbardModelVMC_CL::calc_new_T() const
{
  Eigen::VectorXfp T_new( lat->L );

  for ( cl_uint i = 0; i < lat->L; ++i ) {
    cl_fptype sum = 0.f;
    for ( cl_uint j = 0; j < lat->L; ++j ) {
      sum += v( i, j ) * static_cast<cl_fptype>(
               ( econf.get_site_occ( j ) + econf.get_site_occ( j + lat->L ) ) );
    }
    T_new( i ) = exp( sum );
  }

  return T_new;
}



Eigen::VectorXfp HubbardModelVMC_CL::calc_qupdated_T( const ElectronHop& hop ) const
{
  Eigen::VectorXfp T_prime( lat->L );

  for ( cl_uint i = 0; i < lat->L; ++i ) {
    T_prime( i ) = T( i ) * v.exp( i, lat->get_spinup_site( hop.l ) )
                   / v.exp( i, lat->get_spinup_site( hop.k_pos ) );
  }

  return T_prime;
}



cl_fptype HubbardModelVMC_CL::E_l()
{
  // enqueue transfer of W from the device to the host
  clQ_Wbu.enqueueReadBuffer(
    *devWbu_active, CL_FALSE,
    0, Wbu.size() * sizeof( cl_fptype ), Wbu.data()
  );
  if ( M.ssym == true ) {
    clQ_Wd.enqueueReadBuffer(
      *devWd_active, CL_FALSE,
      0, Wd.size() * sizeof( cl_fptype ), Wd.data()
    );
  }

  // buffer vector for X nearest neighbors
  // (in order to avoid allocating new ones all the time)
  std::vector<cl_uint> k_pos_Xnn;

  // calculate expectation value of the T part of H
  cl_fptype E_l_kin = 0.f;

  // wait for the transfer of W_bu from the device to finish
  clQ_Wbu.finish();

  for ( cl_uint k = 0; k < econf.N(); ++k ) {

    if ( M.ssym == true && k == econf.N() / 2 ) {
      // wait for W_d transfer to finish
      clQ_Wd.finish();
    }

    const cl_uint k_pos = econf.get_electron_pos( k );
    assert( econf.get_site_occ( k_pos ) == ELECTRON_OCCUPATION_FULL );

    for ( cl_uint X = 1; X <= t.size(); ++X ) {
      if ( t[X - 1] == 0.f ) {
        continue;
      }

      cl_fptype sum_Xnn = 0.f;
      lat->get_Xnn( k_pos, X, &k_pos_Xnn );
      for ( auto l_it = k_pos_Xnn.begin(); l_it != k_pos_Xnn.end(); ++l_it ) {
        if ( econf.get_site_occ( *l_it ) == ELECTRON_OCCUPATION_EMPTY ) {
          const cl_fptype R_j = T( lat->get_spinup_site( *l_it ) )
                                / T( lat->get_spinup_site( k_pos ) )
                                * v.exp_onsite() / v.exp( *l_it, k_pos );
          if ( M.ssym == true && k >= econf.N() / 2 ) {
            sum_Xnn += R_j * Wd( *l_it - lat->L, k - econf.N() / 2 );
          } else {
            sum_Xnn += R_j * Wbu( *l_it, k );
          }
        }
      }
      E_l_kin -= t[X - 1] * sum_Xnn;

    }
  }

  const cl_fptype E_l_result =
    ( E_l_kin + U * econf.get_num_dblocc() ) /
    static_cast<cl_fptype>( lat->L );

#if VERBOSE >= 1
  cout << "HubbardModelVMC_CL::E_l() = " << E_l_result << endl;
#endif

  return E_l_result;
}



cl_ulong HubbardModelVMC_CL::mctime() const
{
  return completed_mcsteps;
}



FPDevStat HubbardModelVMC_CL::get_W_devstat() const
{
  return W_devstat;
}
FPDevStat HubbardModelVMC_CL::get_T_devstat() const
{
  return T_devstat;
}
