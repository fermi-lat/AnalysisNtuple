/** @file PointingInfo.cxx
@brief declaration and definition of the class PointingInfo


$Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/PointingInfo.cxx,v 1.13 2007/09/13 22:11:38 burnett Exp $

*/
class MsgStream;
#include "Event/MonteCarlo/Exposure.h"
#include "astro/PointingInfo.h"
#include "AnalysisNtuple/PointingInfo.h"

#include "astro/SkyDir.h"
#include "astro/JulianDate.h"
#include "astro/EarthCoordinate.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

using namespace AnalysisNtuple;

void PointingInfo::set(double time, const astro::PointingInfo& info)
{
    start = time;

    CLHEP::Hep3Vector pos_km = info.position();
    CLHEP::Hep3Vector location = 1.e3* pos_km; 
    
    // cartesian location of the LAT (in m)
    sc_position[0] = location.x();
    sc_position[1] = location.y();
    sc_position[2] = location.z(); 

    ra_zenith = info.zenith().ra();
    dec_zenith= info.zenith().dec();

    ra_scx =    info.xAxis().ra();
    dec_scx =   info.xAxis().dec();

    ra_scz =    info.zAxis().ra();
    dec_scz =   info.zAxis().dec();

    const astro::EarthCoordinate& earth(info.earthCoord());

    lat_geo = earth.latitude(); 
    lon_geo = earth.longitude(); 
    
    rad_geo = earth.altitude(); 
    L=earth.L();
    B=earth.B();
    lat_mag = earth.geolat();
    in_saa= earth.insideSAA()? 1:0;
    zenith_scz = 180/M_PI* info.zenith().difference(info.zAxis());
}
//------------------------------------------------------------------------
Event::Exposure* PointingInfo::forTDS()const
{

    Event::Exposure* entry = new Event::Exposure;
    entry->init(start,lat_geo,lon_geo,rad_geo,
        sc_position[0],sc_position[1],sc_position[2],ra_scx,dec_scx,ra_scz,dec_scz);
    return entry;

}
//------------------------------------------------------------------------
void PointingInfo::finish(double stop_time, double live)
{
    stop = stop_time;
    livetime = live;

}
//------------------------------------------------------------------------
void PointingInfo::setFT2Tuple(INTupleWriterSvc* tuple, const std::string& tname)
{
    if( tuple==0 ) return;
    tuple->addItem(tname, "start",     &start);
    tuple->addItem(tname, "stop",      &stop);
    tuple->addItem(tname, "sc_position[3]", sc_position);
    tuple->addItem(tname, "lat_geo",   &lat_geo);
    tuple->addItem(tname, "lat_mag",   &lat_mag);
    tuple->addItem(tname, "lon_geo",   &lon_geo);
    tuple->addItem(tname, "rad_geo",   &rad_geo);
    tuple->addItem(tname, "ra_zenith", &ra_zenith);
    tuple->addItem(tname, "dec_zenith",&dec_zenith);
    tuple->addItem(tname, "B_McIlwain",&B);
    tuple->addItem(tname, "L_McIlwain",&L);
    tuple->addItem(tname, "ra_scz",    &ra_scz);
    tuple->addItem(tname, "dec_scz",   &dec_scz);
    tuple->addItem(tname, "ra_scx",    &ra_scx);
    tuple->addItem(tname, "dec_scx",   &dec_scx);
    tuple->addItem(tname, "in_saa",   &in_saa);
    tuple->addItem(tname, "livetime",  &livetime);
}
//------------------------------------------------------------------------
void PointingInfo::setPtTuple(INTupleWriterSvc* tuple, const std::string& tname)
{

       /** @page anatup_vars 
     These items are added to the merit tuple  to give the current instrument orientation 


    - PtTime (s) Current time, same as the elapsed time
    - PtLat,PtLon (deg) lattitude and longitude
    - PtAlt (km) altitude
    - PtMagLat magnetic latitude
    - PtPos[3] (m) current orbit position
    - PtRax,PtDecx (deg) equatorial coordinates for orientation of S/C x-axis
    - PtRaz,PtDecz (deg) equatorial coordinates for orientation of S/C z-axis 
    - PtSCzenith  (deg) current angle between zenith and S/C z-axis
    */

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
