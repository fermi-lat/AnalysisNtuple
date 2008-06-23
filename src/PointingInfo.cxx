/** @file PointingInfo.cxx
@brief declaration and definition of the class PointingInfo


$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/PointingInfo.cxx,v 1.3 2008/03/08 18:23:44 lsrea Exp $

*/
class MsgStream;
#include "Event/MonteCarlo/Exposure.h"
#include "astro/PointingInfo.h"
#include "AnalysisNtuple/PointingInfo.h"

#include "astro/SkyDir.h"
#include "astro/GPS.h"
#include "astro/JulianDate.h"
#include "astro/EarthCoordinate.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

using namespace AnalysisNtuple;

void PointingInfo::execute( const astro::GPS& gps)
{
    start = gps.time();
    CLHEP::Hep3Vector pos_km = gps.position();
    CLHEP::Hep3Vector location = 1.e3* pos_km; 
    
    // cartesian location of the LAT (in m)
    sc_position[0] = location.x();
    sc_position[1] = location.y();
    sc_position[2] = location.z(); 

    ra_zenith = gps.zenithDir().ra();
    dec_zenith= gps.zenithDir().dec();

    ra_scx =    gps.xAxisDir().ra();
    dec_scx =   gps.xAxisDir().dec();

    ra_scz =    gps.zAxisDir().ra();
    dec_scz =   gps.zAxisDir().dec();

    const astro::EarthCoordinate& earth(gps.earthpos());

    lat_geo = earth.latitude(); 
    lon_geo = earth.longitude(); 
    
    rad_geo = earth.altitude(); 
    L=earth.L();
    B=earth.B();
    lat_mag = earth.geolat();
    in_saa= earth.insideSAA()? 1:0;
    zenith_scz = 180/M_PI* gps.zenithDir().difference(gps.zAxisDir());
}


/** @page anatup_vars 
    @section Pt  Pt Variables

    These items are added to the merit tuple  to give the current instrument orientation 

<table>
<tr><th> Variable <th>Type<th> Description
<tr><td> PtTime       <td>D<td> (s) Current time, same as the elapsed time
<tr><td> PtLat,PtLon  <td>F<td> (deg) lattitude and longitude
<tr><td> PtAlt        <td>F<td> (km) altitude
<tr><td> PtMagLat     <td>F<td> magnetic latitude
<tr><td> PtPos[3]     <td>F<td> (m) current orbit position
<tr><td> PtRax,PtDecx <td>F<td> (deg) equatorial coordinates for orientation of S/C x-axis
<tr><td> PtRaz,PtDecz <td>F<td> (deg) equatorial coordinates for orientation of S/C z-axis 
<tr><td> PtSCzenith   <td>F<td> (deg) current angle between zenith and S/C z-axis
</table>
*/

//------------------------------------------------------------------------
void PointingInfo::setPtTuple(INTupleWriterSvc* tuple, const std::string& tname)
{
    if( tuple==0 ) return;
    tuple->addItem(tname, "PtTime",   &start);
    tuple->addItem(tname, "PtPos[3]", sc_position);
    tuple->addItem(tname, "PtLat",    &lat_geo);
    tuple->addItem(tname, "PtLon",    &lon_geo);
    tuple->addItem(tname, "PtMagLat", &lat_mag);
    tuple->addItem(tname, "PtAlt",    &rad_geo);
    tuple->addItem(tname, "PtRaz",    &ra_scz);
    tuple->addItem(tname, "PtDecz",   &dec_scz);
    tuple->addItem(tname, "PtRax",    &ra_scx);
    tuple->addItem(tname, "PtDecx",   &dec_scx);
    tuple->addItem(tname, "PtSCzenith", &zenith_scz);
    tuple->addItem(tname, "PtMcIlwainB", &B);
    tuple->addItem(tname, "PtMcIlwainL", &L);

}
