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

#include "mccresults.hpp"

#include <iostream>
#include <fstream>

using namespace std;
namespace fs = boost::filesystem;


void MCCResults::write_to_files( const fs::path& dir ) const
{
  {
    ofstream success_file( ( dir / "sim_success.txt" ).string() );
    success_file << success << endl;
  }

  if ( E ) {
    ofstream E_file( ( dir / "sim_res_E.txt" ).string() );
    E_file << E->mean << " " << E->sigma << endl;
  }

  if ( Deltak ) {
    ofstream Dk_file( ( dir / "sim_res_Dk.txt" ).string() );
    Dk_file << Deltak->transpose() << endl;
  }

  if ( Deltak_Deltakprime ) {
    ofstream DkDkp_file( ( dir / "sim_res_DkDkp.txt" ).string() );
    DkDkp_file << Deltak_Deltakprime.get() << endl;
  }

  if ( Deltak_E ) {
    ofstream DkE_file( ( dir / "sim_res_DkE.txt" ).string() );
    DkE_file << Deltak_E->transpose() << endl;
  }

  if ( dblocc ) {
    ofstream dblocc_file( ( dir / "sim_res_dblocc.txt" ).string() );
    dblocc_file << dblocc->mean << " " << dblocc->sigma << endl;
  }

  if ( nncorr ) {
    ofstream nncorr_file( ( dir / "sim_res_nncorr.txt" ).string() );
    nncorr_file << nncorr.get() << endl;
  }
}


std::ostream& operator<<( std::ostream& out, const MCCResults& res )
{
  out << endl;

  if ( res.success ) {
    out << "Monte Carlo cycle successful!" << endl;
  } else {
    out << "Monte Carlo cycle NOT successful!" << endl;
  }

  if ( res.E ) {
    out << endl
        << "      E = " << res.E->mean << endl
        << "sigma_E = " << res.E->sigma << endl;
  }

  if ( res.Deltak ) {
    out << endl
        << "Delta_k = " << endl
        << res.Deltak->transpose() << endl;
  }

  if ( res.Deltak_Deltakprime ) {
    out << endl
        << "DkDkp = " << endl
        << res.Deltak_Deltakprime.get() << endl;
  }

  if ( res.Deltak_E ) {
    out << endl
        << "Dk_E = " << endl
        << res.Deltak_E->transpose() << endl;
  }

  if ( res.dblocc ) {
    out << endl
        << "      dblocc = " << res.dblocc->mean << endl
        << "sigma_dblocc = " << res.dblocc->sigma << endl;
  }

  if ( res.nncorr ) {
    out << endl
        << "nncorr = " << endl
        << res.nncorr.get() << endl;
  }

  return out;
}
