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

#include "utils.hpp"
using namespace std;
namespace fs = boost::filesystem;


void ostream_setup( ostream& stream )
{
  // stream << setiosflags( ios::scientific );
  stream.setf( ios::showpos );
  stream.precision( numeric_limits<cl_fptype>::digits10 + 1 );
}


fs::path get_hVMC_dir()
{
#if defined(OS_LINUX)
  char buff[1024];
  ssize_t len = readlink( "/proc/self/exe", buff, sizeof( buff ) - 1 );
  if ( len != -1 ) {
    buff[len] = '\0';
    return fs::path( buff ).branch_path();
  } else {
    throw runtime_error( "readlink(/proc/self/exe) failed" );
  }
#elif defined(OS_WINDOWS)
  char buff[MAX_PATH];
  GetModuleFileName( NULL, buff, MAX_PATH );
  return fs::path( buff ).branch_path();
#else
#error "Platform not supported: no way to determine executable directory"
#endif
}


string get_file_contents( const fs::path& file )
{
  ifstream in( file.string(), ios::in | ios::binary );
  if ( in ) {
    string contents;
    in.seekg( 0, ios::end );
    contents.resize( in.tellg() );
    in.seekg( 0, ios::beg );
    in.read( &contents[0], contents.size() );
    in.close();
    return contents;
  } else {
    throw runtime_error( "error while opening the input file" );
  }
}


cl_uint uintsqrt( cl_uint n )
{
  return static_cast<cl_uint>(
           floor( sqrt( static_cast<cl_double>( n ) ) + 0.5 )
         );
}

bool is_perfect_square( cl_uint n )
{
  cl_uint t = uintsqrt( n );
  return t * t == n;
}
