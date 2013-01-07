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
  std::mt19937 rng_init,
  Lattice* const lat_init,
  const SingleParticleOrbitals& M_init,
  const Jastrow& v_init,
  cl_uint N_init,
  const std::vector<cl_fptype>& t_init, cl_fptype U_init,
  cl_fptype W_deviation_target_init,
  cl_uint updates_until_W_recalc_init,
  cl_fptype T_deviation_target_init,
  cl_uint updates_until_T_recalc_init )
  : clCtx( clCtx_init ),
    rng( rng_init ),
    lat( lat_init ), M( M_init ), v( v_init ),
    t( t_init ), U( U_init ),
    econf( ElectronConfiguration( lat, N_init, &rng ) ),
    Wbu(
      M_init.ssym ?
      Eigen::MatrixXfp( lat_init->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat_init->L, N_init )
    ),
    Wd(
      M_init.ssym ?
      Eigen::MatrixXfp( lat_init->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    Wbu_fromdev(
      M_init.ssym ?
      Eigen::MatrixXfp( lat_init->L, N_init / 2 ) :
      Eigen::MatrixXfp( 2 * lat_init->L, N_init )
    ),
    Wd_fromdev(
      M_init.ssym ?
      Eigen::MatrixXfp( lat_init->L, N_init / 2 ) :
      Eigen::MatrixXfp( 0 , 0 )
    ),
    devWbu_active( nullptr ), devWbu_inactive( nullptr ),
    devWd_active( nullptr ), devWd_inactive( nullptr ),
    T( Eigen::VectorXfp( lat_init->L ) ),
    T_fromdev( Eigen::VectorXfp( lat_init->L ) ),
    devT_active( nullptr ), devT_inactive( nullptr ),
    E_l_elbuf( Eigen::VectorXfp( N_init ) ),
    completed_mcsteps( 0 ),
    updates_until_W_recalc( updates_until_W_recalc_init ),
    updates_until_T_recalc( updates_until_T_recalc_init ),
    updates_since_W_recalc( 0 ), updates_since_T_recalc( 0 ),
    W_devstat( FPDevStat( W_deviation_target_init ) ),
    T_devstat( FPDevStat( T_deviation_target_init ) )
{
  // ----- setup of all OpenCL objects -----

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
  stringstream compileropts;
//  compileropts << "-g ";
  compileropts <<  "-D LATTICE_TYPE=" << static_cast<cl_uint>( lat->type );
  compileropts << " -D LATTICE_NUM_SITES=" << lat->L;
  if ( lat->type == LATTICE_2DSQUARE ) {
    compileropts << " -D LATTICE_SIZE=" << uintsqrt( lat->L );
  }
  compileropts << " -D NUM_ELECTRONS=" << econf.N();
  if ( M.ssym == true ) {
    compileropts << " -D M_IS_SPIN_SYMMETRIC";
  }
#ifdef USE_FP_DBLPREC
  compileropts << " -D USE_FP_DBLPREC_OPENCL";
#endif
#if VERBOSE >= 1
  cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : building cl::Program with "
       << compileropts.str() << endl;
#endif
  try {
    clPrg.build( used_devices, compileropts.str().c_str() );
  } catch ( cl::Error& e ) {
    cerr << "Build failed! " << e.what() << " (" << e.err() << ")" << endl;
    cerr << "Build log:" << endl;
    cerr << clPrg.getBuildInfo<CL_PROGRAM_BUILD_LOG>( used_devices[0] ) << endl;
    throw e;
  }

  // set up Kernel objects
  clK_hop           = cl::Kernel( clPrg, "hop" );
  clK_update_devWbu = cl::Kernel( clPrg, "update_Wbu" );
  if ( M.ssym == true ) {
    clK_update_devWd  = cl::Kernel( clPrg, "update_Wd" );
  }
  clK_update_devT   = cl::Kernel( clPrg, "update_T" );
  clK_calc_E_l      = cl::Kernel( clPrg, "calc_E_l" );

  // setup the command queue
  clQ = cl::CommandQueue(
          clCtx, used_devices[0]
        );

  // set up and transfer the Jastrow
  vector<cl_fptype> expv = v.get_reduced_raw_jastrow();
  devexpv = cl::Buffer(
              clCtx, CL_MEM_READ_ONLY,
              expv.size() * sizeof( cl_fptype )
            );
  clQ.enqueueWriteBuffer(
    devexpv, CL_FALSE,
    0, expv.size() * sizeof( cl_fptype ), expv.data()
  );

  // set up and transfer hopping parameters + U buffer
  devUt = cl::Buffer(
            clCtx, CL_MEM_READ_ONLY,
            4 * sizeof( cl_fptype )
          );

  cl_fptype Ut[4];
  Ut[0] = U;
  Ut[1] = t[0];
  Ut[2] = t[1];
  Ut[3] = t[2];

  clQ.enqueueWriteBuffer(
    devUt, CL_FALSE,
    0, 4 * sizeof( cl_fptype ), Ut
  );

  // set up the device buffer for the current electron hop
  devehop = cl::Buffer( clCtx, CL_MEM_READ_WRITE, 128 );

  // set up electronic configuration buffers on the device
  deveconf_site_occ = cl::Buffer(
                        clCtx, CL_MEM_READ_WRITE,
                        2 * lat->L * sizeof( cl_uint )
                      );
  deveconf_electron_pos = cl::Buffer(
                            clCtx, CL_MEM_READ_WRITE,
                            lat->L * sizeof( cl_uint )
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

  // set up T buffers on the device
  devT_1 = cl::Buffer(
             clCtx, CL_MEM_READ_WRITE,
             T.size() * sizeof( cl_fptype )
           );
  devT_2 = cl::Buffer(
             clCtx, CL_MEM_READ_WRITE,
             T.size() * sizeof( cl_fptype )
           );
  devT_active   = &devT_1;
  devT_inactive = &devT_2;

  // set up the E_l output buffer for the individual electrons
  devE_l_elbuf = cl::Buffer(
                   clCtx, CL_MEM_WRITE_ONLY,
                   E_l_elbuf.size() * sizeof( cl_fptype )
                 );

  // ----- initial state generation -----

  bool enough_overlap;

  do {

    econf.distribute_random();
    enough_overlap = true; // assume until we are proven wrong

#if VERBOSE >=1
    cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : checking newly generated "
         << "state for enough overlap" << endl;
#endif

    cl_fptype Wbu_avg, Wd_avg;

    if ( M.ssym == true ) {

      // check spin up part
      Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Du().transpose() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : spin up part has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      Wbu.noalias() = lu_decomp.solve( M.orbitals.transpose() ).transpose();
      Wbu_avg = Wbu.squaredNorm() / static_cast<cl_fptype>( Wbu.size() );
      enough_overlap &= Wbu_avg < 100.f;
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : spin up part has too "
             << "little overlap with the determinantal wavefunction, "
             << "inverse overlap measure is: " << Wbu_avg << endl;
#endif
        continue;
      }

      // check spin down part
      lu_decomp.compute( calc_Dd().transpose() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : spin down part has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      Wd.noalias() = lu_decomp.solve( M.orbitals.transpose() ).transpose();
      Wd_avg = Wd.squaredNorm()  / static_cast<cl_fptype>( Wd.size() );
      enough_overlap &= Wd_avg < 100.f;
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : spin down part has too "
             << "little overlap with the determinantal wavefunction, "
             << "inverse overlap measure is: " << Wd_avg << endl;
#endif
        continue;
      }

    } else {

      // check whole determinantal part
      Eigen::FullPivLU<Eigen::MatrixXfp> lu_decomp( calc_Db() );
      enough_overlap &= lu_decomp.isInvertible();
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : state has no "
             << "overlap with the determinantal wavefunction" << endl;
#endif
        continue;
      }

      Wbu.noalias() = lu_decomp.solve( M.orbitals.transpose() ).transpose();
      Wbu_avg = Wbu.squaredNorm() / static_cast<cl_fptype>( Wbu.size() );
      enough_overlap &= Wbu_avg < 50.f;
      if ( !enough_overlap ) {
#if VERBOSE >= 1
        cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : state has too "
             << "little overlap with the determinantal wavefunction, "
             << "inverse overlap measure is: " << Wbu_avg << endl;
#endif
        continue;
      }

    }

    // check Jastrow part
    calc_new_T();
    cl_fptype T_avg = T.squaredNorm() / static_cast<cl_fptype>( T.size() );
    enough_overlap &= T_avg < 100.f;
    if ( !enough_overlap ) {
#if VERBOSE >= 1
      cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : Jastrow ratios "
           << "are to small, inverse measure is: " << T_avg << endl;
#endif
      continue;
    }

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : state has sufficient "
         << "overlap! inverse overlap measures are: "
         << Wbu_avg << " " << Wd_avg << " " << T_avg
         << " -> initial state selection completed!" << endl;
#endif

  } while ( !enough_overlap );

#if VERBOSE >= 1
  cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : calculated initial matrix W ="
       << endl << Wbu << endl;
  if ( M.ssym == true ) {
    cout << "----->" << endl << Wd << endl;
  }
  cout << "HubbardModelVMC_CL::HubbardModelVMC_CL() : calculated initial vector T ="
       << endl << T << endl;
#endif


  // ----- transfer of initial state to the OpenCL device -----

  // enqueue transfer of Wbu and Wd to the device
  clQ.enqueueWriteBuffer(
    *devWbu_active, CL_FALSE,
    0, Wbu.size() * sizeof( cl_fptype ), Wbu.data()
  );
  if ( M.ssym == true ) {
    clQ.enqueueWriteBuffer(
      *devWd_active, CL_FALSE,
      0, Wd.size() * sizeof( cl_fptype ), Wd.data()
    );
  }

  // enqueue transfer of T to the device
  clQ.enqueueWriteBuffer(
    *devT_active, CL_FALSE,
    0, T.size() * sizeof( cl_fptype ), T.data()
  );

  // copy the current electron configuration to the raw data buffers
  vector<cl_uint> econf_site_occ     = econf.get_site_occ_raw();
  vector<cl_uint> econf_electron_pos = econf.get_electron_pos_raw();

  // enqueue transfer of the electron position buffers to the device
  clQ.enqueueWriteBuffer(
    deveconf_site_occ, CL_FALSE,
    0, econf_site_occ.size() * sizeof( cl_uint ), econf_site_occ.data()
  );
  clQ.enqueueWriteBuffer(
    deveconf_electron_pos, CL_FALSE,
    0, econf_electron_pos.size() * sizeof( cl_uint ), econf_electron_pos.data()
  );

  // wait for all transfers to the device to complete
  clQ.finish();
}



HubbardModelVMC_CL::~HubbardModelVMC_CL()
{
  delete lat;
}



void HubbardModelVMC_CL::sync_econf_down()
{
  vector<cl_uint> electron_pos_raw( econf.N() );

  // download the electron position buffers from the device (blocking!)
  clQ.enqueueReadBuffer(
    deveconf_electron_pos, CL_TRUE,
    0, electron_pos_raw.size() * sizeof( cl_uint ), electron_pos_raw.data()
  );

  // use the raw position vector to make a new ElectronConfiguration object
  econf.init_from_raw_elpos( electron_pos_raw );
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



void HubbardModelVMC_CL::metstep()
{
  // generate the random numbers
  cl_uint hopping_elid
    = uniform_int_distribution<cl_uint>( 0, econf.N() - 1 )( rng );
  cl_uint hopto_nid
    = uniform_int_distribution<cl_uint>( 0, lat->coordnum - 1 )( rng );
  cl_fptype hop_selector
    = uniform_real_distribution<cl_fptype>( 0, 1 )( rng ); 

  // enqueue a new hop
  clK_hop.setArg( 0, devehop );
  clK_hop.setArg( 1, deveconf_electron_pos );
  clK_hop.setArg( 2, deveconf_site_occ );
  clK_hop.setArg( 3, *devWbu_active );
  clK_hop.setArg( 4, *devWd_active );
  clK_hop.setArg( 5, *devT_active );
  clK_hop.setArg( 6, devexpv );
  clK_hop.setArg( 7, hopping_elid );
  clK_hop.setArg( 8, hopto_nid );
  clK_hop.setArg( 9, hop_selector );
  clQ.enqueueNDRangeKernel( clK_hop,
    cl::NullRange, cl::NDRange( 1 ), cl::NDRange( 1 )
  );

  // TODO: remove (Robert Rueger, 2013-01-06 15:22)
  clQ.finish();

  // update W and T
  perform_W_update();
  perform_T_update();
}



void HubbardModelVMC_CL::perform_W_update()
{
  if ( updates_since_W_recalc >= updates_until_W_recalc ) {

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_W_update() : recalculating W!" << endl;
#endif

    // download electron configuration from the device
    sync_econf_down();

    // enqueue update of W on the device
    // (we want to update it to the point of the recalculated W)
    // (the updated W will be in the _inactive_ buffer!)
    calc_qupdated_W();

    // enqueue transfer of the updated W from the device to the host
    // (waits for the update to finish first)
    clQ.enqueueReadBuffer(
      *devWbu_inactive, CL_FALSE,
      0, Wbu_fromdev.size() * sizeof( cl_fptype ), Wbu_fromdev.data()
    );
    if ( M.ssym == true ) {
      clQ.enqueueReadBuffer(
        *devWd_inactive, CL_FALSE,
        0, Wd_fromdev.size() * sizeof( cl_fptype ), Wd_fromdev.data()
      );
    }

    // recalculate W on the host
    calc_new_W();

    // wait for transfer of W from the device to finish
    clQ.finish();

    // enqueue upload of the recalculated W to the device
    // (it will be put in the _active_ buffer)
    clQ.enqueueWriteBuffer(
      *devWbu_active, CL_FALSE,
      0, Wbu.size() * sizeof( cl_fptype ), Wbu.data()
    );
    if ( M.ssym == true ) {
      clQ.enqueueWriteBuffer(
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

    if ( !( dev < W_devstat.target ) ) {
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

    calc_qupdated_W();

    swap( devWbu_active, devWbu_inactive );
    if ( M.ssym == true ) {
      swap( devWd_active, devWd_inactive );
    }

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

    Wbu.noalias()
      = calc_Du().transpose()
        .partialPivLu().solve( M.orbitals.transpose() )
        .transpose();
    Wd.noalias()
      = calc_Dd().transpose()
        .partialPivLu().solve( M.orbitals.transpose() )
        .transpose();

  } else {

    Wd.noalias()
      = calc_Db().transpose()
        .partialPivLu().solve( M.orbitals.transpose() )
        .transpose();

  }
}



void HubbardModelVMC_CL::calc_qupdated_W()
{
  // enqueue update of Wbu
  clK_update_devWbu.setArg( 0, devehop );
  clK_update_devWbu.setArg( 1, *devWbu_active );
  clK_update_devWbu.setArg( 2, *devWbu_inactive );
  clQ.enqueueNDRangeKernel(
    clK_update_devWbu,
    cl::NullRange, cl::NDRange( Wbu.size() ), cl::NullRange
  );

  if ( M.ssym == true ) {
    // enqueue update of Wd
    clK_update_devWd.setArg( 0, devehop );
    clK_update_devWd.setArg( 1, *devWd_active );
    clK_update_devWd.setArg( 2, *devWd_inactive );
    clQ.enqueueNDRangeKernel(
      clK_update_devWd,
      cl::NullRange, cl::NDRange( Wd.size() ), cl::NullRange
    );
  }
}



void HubbardModelVMC_CL::perform_T_update()
{
  if ( updates_since_T_recalc >= updates_until_T_recalc ) {

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_T_update() : recalculating T!" << endl;
#endif

    // download electron configuration from the device
    sync_econf_down();
    
    // enqueue update of T on the device
    // (we want to update it to the point of the recalculated T)
    // (the updated T will be in the _inactive_ buffer!)
    calc_qupdated_T();

    // enqueue transfer of the updated T from the device to the host
    // (waits for the update to finish first)
    clQ.enqueueReadBuffer(
      *devT_inactive, CL_FALSE,
      0, T_fromdev.size() * sizeof( cl_fptype ), T_fromdev.data()
    );

    // recalculate T on the host
    calc_new_T();

    // wait for transfer of T from the device to finish
    clQ.finish();

    // enqueue upload of the recalculated T to the device
    // (it will be put in the _active_ buffer)
    clQ.enqueueWriteBuffer(
      *devT_active, CL_FALSE,
      0, T.size() * sizeof( cl_fptype ), T.data()
    );

    // compare the recalculated T and the updated T on the host
    cl_fptype dev = calc_deviation( T, T_fromdev );
    T_devstat.add( dev );

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_T_update() : recalculated T "
         << "with deviation = " << dev << endl;

    if ( !( dev < T_devstat.target ) ) {
      cout << "HubbardModelVMC_CL::perform_T_update() : deviation goal for matrix "
           << "T not met!" << endl
           << "HubbardModelVMC_CL::perform_T_update() : approximate T =" << endl
           << T_fromdev << endl;
      cout << "HubbardModelVMC_CL::perform_T_update() : exact T =" << endl
           << T << endl;
    }
#endif

    assert( dev < T_devstat.target );

    updates_since_T_recalc = 0;

  } else {

#if VERBOSE >= 1
    cout << "HubbardModelVMC_CL::perform_T_update() : "
         << "performing a quick update of T!" << endl;
#endif

    calc_qupdated_T();

    swap( devT_active, devT_inactive );

    ++updates_since_T_recalc;
  }
}



void HubbardModelVMC_CL::calc_new_T()
{
  for ( cl_uint i = 0; i < lat->L; ++i ) {
    cl_fptype sum = 0.f;
    for ( cl_uint j = 0; j < lat->L; ++j ) {
      sum += v( i, j ) * static_cast<cl_fptype>(
               ( econf.get_site_occ( j ) + econf.get_site_occ( j + lat->L ) ) );
    }
    T( i ) = exp( sum );
  }
}



void HubbardModelVMC_CL::calc_qupdated_T()
{
  clK_update_devT.setArg( 0, devehop );
  clK_update_devT.setArg( 1, *devT_active );
  clK_update_devT.setArg( 2, *devT_inactive );
  clK_update_devT.setArg( 3, devexpv );
  clQ.enqueueNDRangeKernel(
    clK_update_devT,
    cl::NullRange, cl::NDRange( lat->L ), cl::NullRange
  );
}



cl_fptype HubbardModelVMC_CL::E_l()
{
/*
  // enqueue calculation of E_l for the individual electrons
  clK_calc_E_l.setArg( 0, devE_l_elbuf );
  clK_calc_E_l.setArg( 1, deveconf_electron_pos );
  clK_calc_E_l.setArg( 2, deveconf_site_occ );
  clK_calc_E_l.setArg( 3, *devWbu_active );
  clK_calc_E_l.setArg( 4, *devWd_active );
  clK_calc_E_l.setArg( 5, *devT_active );
  clK_calc_E_l.setArg( 6, devexpv );
  clK_calc_E_l.setArg( 7, devUt );
  clQ.enqueueNDRangeKernel( clK_calc_E_l,
    cl::NullRange, cl::NDRange( econf.N() ), cl::NullRange
  );

  // transfer results to the host (blocking!)
  clQ.enqueueReadBuffer(
    devE_l_elbuf, CL_TRUE,
    0, E_l_elbuf.size() * sizeof(cl_fptype), E_l_elbuf.data()
  );

  // return the summed up results of the individual electrons
  cout << E_l_elbuf.sum() / static_cast<cl_fptype>( lat->L ) << endl;
  return E_l_elbuf.sum() / static_cast<cl_fptype>( lat->L );
*/

  // download electron configuration from the device
  sync_econf_down();

  // enqueue transfer of W from the device to the host
  clQ.enqueueReadBuffer(
    *devWbu_inactive, CL_FALSE,
    0, Wbu.size() * sizeof( cl_fptype ), Wbu.data()
  );
  if ( M.ssym == true ) {
    clQ.enqueueReadBuffer(
      *devWd_inactive, CL_FALSE,
      0, Wd.size() * sizeof( cl_fptype ), Wd.data()
    );
  }

  // enqueue transfer of from the device to the host
  clQ.enqueueReadBuffer(
    *devT_inactive, CL_FALSE,
    0, T.size() * sizeof( cl_fptype ), T.data()
  );

  // wait for all transfers to complete
  clQ.finish();

  // buffer vector for X nearest neighbors
  // (in order to avoid allocating new ones all the time)
  std::vector<cl_uint> k_pos_Xnn;

  // calculate expectation value of the T part of H
  cl_fptype E_l_kin = 0.f;

  for ( cl_uint k = 0; k < econf.N(); ++k ) {

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
