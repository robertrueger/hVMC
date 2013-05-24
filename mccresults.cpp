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

#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "serialization_eigen.hpp"

using namespace std;
namespace fs = boost::filesystem;
namespace ar  = boost::archive;


void MCCResults::write_to_files( const fs::path& dir ) const
{
  if ( E ) {
    ofstream E_txtfile( ( dir / "sim_res_E.txt" ).string() );
    E_txtfile << E->mean << " " << E->sigma << endl;

    ofstream E_datfile( ( dir / "sim_res_E.dat" ).string() );
    ar::text_oarchive E_archive( E_datfile );
    E_archive << E.get();
  }

  if ( Deltak ) {
    ofstream Dk_txtfile( ( dir / "sim_res_Dk.txt" ).string() );
    Dk_txtfile << Deltak->transpose() << endl;

    ofstream Dk_datfile( ( dir / "sim_res_Dk.dat" ).string() );
    ar::text_oarchive Dk_archive( Dk_datfile );
    Dk_archive << Deltak.get();
  }

  if ( Deltak_Deltakprime ) {
    ofstream DkDkp_file( ( dir / "sim_res_DkDkp.txt" ).string() );
    DkDkp_file << Deltak_Deltakprime.get() << endl;

    ofstream DkDkp_datfile( ( dir / "sim_res_DkDkp.dat" ).string() );
    ar::text_oarchive DkDkp_archive( DkDkp_datfile );
    DkDkp_archive << Deltak_Deltakprime.get();
  }

  if ( Deltak_E ) {
    ofstream DkE_txtfile( ( dir / "sim_res_DkE.txt" ).string() );
    DkE_txtfile << Deltak_E->transpose() << endl;

    ofstream DkE_datfile( ( dir / "sim_res_DkE.dat" ).string() );
    ar::text_oarchive DkE_archive( DkE_datfile );
    DkE_archive << Deltak_E.get();
  }

  if ( dblocc ) {
    ofstream dblocc_txtfile( ( dir / "sim_res_dblocc.txt" ).string() );
    dblocc_txtfile << dblocc->mean << " " << dblocc->sigma << endl;

    ofstream dblocc_datfile( ( dir / "sim_res_dblocc.dat" ).string() );
    ar::text_oarchive dblocc_archive( dblocc_datfile );
    dblocc_archive << dblocc.get();
  }

  if ( nncorr ) {
    ofstream nncorr_txtfile( ( dir / "sim_res_nncorr.txt" ).string() );
    nncorr_txtfile << nncorr.get() << endl;

    ofstream nncorr_datfile( ( dir / "sim_res_nncorr.dat" ).string() );
    ar::text_oarchive nncorr_archive( nncorr_datfile );
    nncorr_archive << nncorr.get();
  }

  if ( sscorr ) {
    ofstream sscorr_txtfile( ( dir / "sim_res_sscorr.txt" ).string() );
    sscorr_txtfile << sscorr.get() << endl;

    ofstream sscorr_datfile( ( dir / "sim_res_sscorr.dat" ).string() );
    ar::text_oarchive sscorr_archive( sscorr_datfile );
    sscorr_archive << sscorr.get();
  }

  if ( pconfs ) {
    ofstream pconfs_txtfile( ( dir / "sim_res_pconfs.txt" ).string() );
    for ( auto it = pconfs.get().begin(); it != pconfs.get().end(); ++it ) {
      pconfs_txtfile << it->transpose() << endl;
    }

    ofstream pconfs_datfile( ( dir / "sim_res_pconfs.dat" ).string() );
    ar::text_oarchive pconfs_archive( pconfs_datfile );
    pconfs_archive << pconfs.get();
  }
}


std::ostream& operator<<( std::ostream& out, const MCCResults& res )
{
  out << endl;

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

  if ( res.sscorr ) {
    out << endl
        << "sscorr = " << endl
        << res.sscorr.get() << endl;
  }

  if ( res.pconfs ) {
    cout << endl
         << "pconfs = " << endl;
    for ( auto it = res.pconfs.get().begin();
          it != res.pconfs.get().end();
          ++it ) {
      cout << it->transpose() << endl;
    }
  }

  out << endl;

  return out;
}
