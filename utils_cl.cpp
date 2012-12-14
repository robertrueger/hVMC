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

#include "utils_cl.hpp"
using namespace std;


void print_clinfo()
{
  cout << endl;

  vector<cl::Platform> platform_list;
  cl::Platform::get( &platform_list );

  cout << "Discovered " << platform_list.size() << " OpenCL platform";
  if ( platform_list.size() != 1 ) {
    cout << "s";
  }
  cout << endl << endl;

  cl_uint pl_count = 0;
  for ( auto pl = platform_list.begin(); pl != platform_list.end(); ++pl ) {

    cout << "=======> PLATFORM #" << pl_count << endl;

    string pl_name;
    pl->getInfo( CL_PLATFORM_NAME, &pl_name );
    cout << "   Name: " << pl_name << endl;

    string pl_vendor;
    pl->getInfo( CL_PLATFORM_VENDOR, &pl_vendor );
    cout << " Vendor: " << pl_vendor << endl;

    string pl_version;
    pl->getInfo( CL_PLATFORM_VERSION, &pl_version );
    cout << "Version: " << pl_version << endl;

    string pl_profile;
    pl->getInfo( CL_PLATFORM_PROFILE, &pl_profile );
    cout << "Profile: " << pl_profile << endl;

    vector<cl::Device> pl_device_list;
    pl->getDevices( CL_DEVICE_TYPE_ALL, &pl_device_list );
    cout << "Devices: " << pl_device_list.size() << endl;

    cl_uint dev_count = 0;
    for ( auto dev = pl_device_list.begin(); dev != pl_device_list.end(); ++dev ) {

      cout << "---------------> DEVICE #" << dev_count << endl;

      string dev_name;
      dev->getInfo( CL_DEVICE_NAME, &dev_name );
      cout << "           Name: " << dev_name << endl;

      string dev_vendor;
      dev->getInfo( CL_DEVICE_VENDOR, &dev_vendor );
      cout << "         Vendor: " << dev_vendor << endl;

      cl_device_type dev_type;
      dev->getInfo( CL_DEVICE_TYPE, &dev_type );
      cout << "           Type: ";
      if ( dev_type & CL_DEVICE_TYPE_GPU ) {
        cout << "GPU";
      } else if ( dev_type & CL_DEVICE_TYPE_CPU ) {
        cout << "CPU";
      } else {
        cout << "other";
      }
      cout << endl;

      string dev_extensions;
      dev->getInfo( CL_DEVICE_EXTENSIONS, &dev_extensions );
      cout << "    cl_khr_fp64: ";
      if ( dev_extensions.find( "cl_khr_fp64" ) != string::npos ) {
        cout << "supported";
      } else {
        cout << "not supported";
      }
      cout << endl;

      cl_ulong dev_gmemsize;
      dev->getInfo( CL_DEVICE_GLOBAL_MEM_SIZE, &dev_gmemsize );
      cout << "      GMEM size: " << dev_gmemsize << endl;

      cl_bool dev_gmemuh;
      dev->getInfo( CL_DEVICE_HOST_UNIFIED_MEMORY, &dev_gmemuh );
      cout << "GMEM = host mem: ";
      if ( dev_gmemuh == CL_TRUE ) {
        cout << "yes";
      } else {
        cout << "no";
      }
      cout << endl;

      cl_device_mem_cache_type dev_gmemct;
      dev->getInfo( CL_DEVICE_GLOBAL_MEM_CACHE_TYPE, &dev_gmemct );
      cout << "     GMEM cache: ";
      if ( dev_gmemct == CL_NONE ) {
        cout << "none";
      } else if ( dev_gmemct == CL_READ_ONLY_CACHE ) {
        cout << "read only";
      } else { /* ( dev_gmemct == CL_READ_WRITE_CACHE ) */
        cout << "read/write";
      }
      cout << endl;

      cl_ulong dev_gmemcs;
      dev->getInfo( CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, &dev_gmemcs );
      cout << "GMEM cache size: " << dev_gmemcs << endl;

      cl_device_local_mem_type dev_lmemt;
      dev->getInfo( CL_DEVICE_LOCAL_MEM_TYPE, &dev_lmemt );
      cout << "      LMEM type: ";
      if ( dev_lmemt == CL_LOCAL ) {
        cout << "dedicated";
      } else if ( dev_lmemt == CL_GLOBAL ) {
        cout << "global";
      } else if ( dev_lmemt == CL_NONE ) {
        cout << "none";
      } else {
        cout << "unknown";
      }
      cout << endl;

      cl_ulong dev_lmems;
      dev->getInfo( CL_DEVICE_LOCAL_MEM_SIZE, &dev_lmems );
      cout << "      LMEM size: " << dev_lmems << endl;

      cl_ulong dev_cmems;
      dev->getInfo( CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, &dev_cmems );
      cout << "      CMEM size: " << dev_cmems << endl;

      cl_uint dev_compunits;
      dev->getInfo( CL_DEVICE_MAX_COMPUTE_UNITS, &dev_compunits );
      cout << "  Compute units: " << dev_compunits << endl;

      size_t dev_wgs;
      dev->getInfo( CL_DEVICE_MAX_WORK_GROUP_SIZE, &dev_wgs );
      vector<size_t> dev_maxwis;
      dev->getInfo( CL_DEVICE_MAX_WORK_ITEM_SIZES, &dev_maxwis );
      cout << " Max group size: " << dev_wgs << " (per dim: ";
      copy( dev_maxwis.begin(), dev_maxwis.end(),
            ostream_iterator<size_t>( cout, " " ) );
      cout << "\b)" << endl;

      cl_uint dev_vints, dev_vfloats, dev_vdoubles;
      dev->getInfo( CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, &dev_vints );
      dev->getInfo( CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, &dev_vfloats );
      dev->getInfo( CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, &dev_vdoubles );
      cout << "   Vector width: "
           << "int" << dev_vints << " "
           << "float" << dev_vfloats;
      if ( dev_extensions.find( "cl_khr_fp64" ) != string::npos ) {
        cout << " double" << dev_vdoubles << endl;
      }

      ++dev_count;
    }

    ++pl_count;
    cout << endl;
  }
}



cl::Context clcontext_setup( cl_uint pl_id, cl_uint dev_id )
{
  vector<cl::Device> selected_device;

  try {

    vector<cl::Platform> platform_list;
    cl::Platform::get( &platform_list );

    if ( pl_id > platform_list.size() - 1 ) {
      throw runtime_error( "specified platform not found" );
    }

    vector<cl::Device> device_list;
    platform_list.at( pl_id ).getDevices( CL_DEVICE_TYPE_ALL, &device_list );

    if ( dev_id > device_list.size() - 1 ) {
      throw runtime_error( "specified device not found" );
    }

    selected_device.push_back( device_list.at( dev_id ) );

#ifdef USE_FP_DBLPREC
    string dev_extensions;
    selected_devices.at( 0 ).getInfo( CL_DEVICE_EXTENSIONS, &dev_extensions );
    if ( dev_extensions.find( "cl_khr_fp64" ) == string::npos ) {
      throw runtime_error( "selected device does not support double precision" );
    }
#endif

  } catch ( const runtime_error& e ) {
    cerr << "Error while creating OpenCL context: " << e.what() << endl;
    exit( 1 );
  }

  return cl::Context( selected_device );
}
