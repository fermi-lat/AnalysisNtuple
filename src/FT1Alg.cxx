/** @file FT1Alg.cxx
@brief Declaration and implementation of Gaudi algorithm FT1Alg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/FT1Alg.cxx,v 1.4 2007/06/07 17:00:13 lsrea Exp $
*/
// Include files

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "astro/PointingTransform.h"

#include "Event/MonteCarlo/Exposure.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <cassert>
#include <map>

// forward declatation of the worker
class FT1worker;

namespace { // anonymous namespace for file-global
    unsigned int nbOfEvtsInFile(100000);
    std::string treename("MeritTuple");
#include "Item.h"
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class FT1Alg
@brief Extract info from tuple, etc. to add ft1 items to this of another tree
*/
class FT1Alg : public Algorithm {

public:
    FT1Alg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
private:
    /// this guy does the work!
    FT1worker * m_worker;
    //counter
    int m_count;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static const AlgFactory<FT1Alg>  Factory;
const IAlgFactory& FT1AlgFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class FT1worker{ 
public:
    FT1worker();

    void evaluate(const Event::Exposure* exp);

private:
    std::map<std::string, double> getCelestialCoords(const Event::Exposure& exp,
        const CLHEP::Hep3Vector glastDir);

    bool useVertex(){ //TODO: implement
        return false;
    }

    // tuple items expect to find
    TypedItem<unsigned int, 'i'> EvtRun, EvtEventId;

    // these all float or double
    Item EvtLiveTime;
    Item EvtEnergyCorr;
    Item EvtElapsedTime;
    Item VtxXDir, VtxYDir, VtxZDir;
    Item VtxX0, VtxY0, VtxZ0;
    Item TkrNumTracks;
    Item Tkr1XDir, Tkr1YDir, Tkr1ZDir;
    Item Tkr1X0, Tkr1Y0, Tkr1Z0;
    Item Tkr1FirstLayer;
    Item CTBBestEnergy;
    Item CTBBestXDir;
    Item CTBBestYDir;
    Item CTBBestZDir;

    //FT1 entries to create
    unsigned int m_ft1eventid;
    float m_ft1energy;
    float m_ft1theta,m_ft1phi,m_ft1ra,m_ft1dec,m_ft1l,m_ft1b;
    float m_ft1zen,m_ft1azim,m_ft1livetime;
    float m_ft1convpointx,m_ft1convpointy,m_ft1convpointz,m_ft1convlayer;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FT1Alg::FT1Alg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)
{
    declareProperty("TreeName",  treename="MeritTuple");
    declareProperty("NbOfEvtsInFile", nbOfEvtsInFile=100000);

}

StatusCode FT1Alg::initialize()
{
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // Use the Job options service to get the Algorithm's parameters
    setProperties();

    // get a pointer to RootTupleSvc 
    if( (sc = service("RootTupleSvc", rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
        return sc;
    }
    m_worker = new FT1worker();

    m_count = 0;

    return sc;
}

StatusCode FT1Alg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    m_count++;
    //First get the coordinates from the ExposureCol

    const Event::Exposure* exp(0);

    Event::ExposureCol* elist = 0;
    eventSvc()->retrieveObject("/Event/MC/ExposureCol",(DataObject *&)elist);
    if ( elist==0 ) {
#if 1

#else
        if(m_count<6) {
            log << MSG::INFO << "No ExposureCol found" << endreq;
        }
        else if (m_count==6) {
            log << MSG::INFO << "Message suppressed after 5 events" << endreq;
        }
#endif
    }else {
        exp = *(*elist).begin();
    }
    // now have the worker do it
    m_worker->evaluate(exp);

    return sc;
}

StatusCode FT1Alg::finalize()
{
    return StatusCode::SUCCESS;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FT1worker::FT1worker()
// initialize pointers to current items
: EvtRun("EvtRun")
, EvtEventId("EvtEventId")
, EvtElapsedTime("EvtElapsedTime")
, EvtLiveTime("EvtLiveTime")
, EvtEnergyCorr("EvtEnergyCorr")
, TkrNumTracks("TkrNumTracks")
, VtxXDir("VtxXDir")
, VtxYDir("VtxYDir")
, VtxZDir("VtxZDir")
, VtxX0("VtxX0")
, VtxY0("VtxY0")
, VtxZ0("VtxZ0")
, Tkr1XDir("Tkr1XDir")
, Tkr1YDir("Tkr1YDir")
, Tkr1ZDir("Tkr1ZDir")
, Tkr1X0("Tkr1X0")
, Tkr1Y0("Tkr1Y0")
, Tkr1Z0("Tkr1Z0")
, Tkr1FirstLayer("Tkr1FirstLayer")
, CTBBestEnergy("CTBBestEnergy")
, CTBBestXDir("CTBBestXDir")  
, CTBBestYDir("CTBBestYDir")  
, CTBBestZDir("CTBBestZDir")
{
    //now create new items 

    addItem( "FT1EventId",       m_ft1eventid);
    addItem( "FT1Energy",        m_ft1energy);
    addItem( "FT1Theta",         m_ft1theta);
    addItem( "FT1Phi",           m_ft1phi);
    addItem( "FT1Ra",            m_ft1ra);
    addItem( "FT1Dec",           m_ft1dec);
    addItem( "FT1L",             m_ft1l);
    addItem( "FT1B",             m_ft1b);
    addItem( "FT1ZenithTheta",   m_ft1zen);
    addItem( "FT1EarthAzimuth",  m_ft1azim);
    addItem( "FT1ConvPointX",    m_ft1convpointx);
    addItem( "FT1ConvPointY",    m_ft1convpointy);
    addItem( "FT1ConvPointZ",    m_ft1convpointz);
    addItem( "FT1ConvLayer",     m_ft1convlayer);
    addItem( "FT1Livetime",      m_ft1livetime);
}

/** @page anatup_vars 

@section FT1  FT1 Variables

see <a href="http://glast.gsfc.nasa.gov/ssc/dev/fits_def/definitionFT1.html">FT1 definition</a>

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> FT1EventId  
<td>UI<td>  RunNumber*100000 + EventNumber  
<tr><td> FT1Energy   
<td>F<td>  (MeV) estimate for energy  
<tr><td> FT1Theta, FT1Phi 
<td>F<td>  (deg) reconstructed direction with respect to instrument coordinate system      
<tr><td> FT1Ra, FT1Dec 
<td>F<td>  (deg) reconstructed direction in equatorial coordinates       
<tr><td> FT1ZenithTheta, FT1EarthAzimuth 
<td>F<td>  (deg) reconstucted direction with respect to local zenith system
<tr><td> FT1L, FT1B 
<td>F<td>  (deg) galactic longitude and latitude of reconstructed direction
<tr><td> FT1Livetime 
<td>F<td>  (s) Cumulative live time, from start of run, or mission
<tr><td> FT1ConvLayer 
<td>F<td>  Starting layer of the best track found in the tracker 
(Layer 0 is the one closest to the calorimeter.)
<tr><td> FT1ConvPoint[X/Y/Z] 
<td>F<td>  <b>Do not use; no longer filled!</b>
</table> 
*/

void FT1worker::evaluate(const Event::Exposure* exp)
{


    //eventId and Time are always defined 
    m_ft1eventid =  nbOfEvtsInFile *EvtRun.value()  + EvtEventId.value();

    // Give default "guard" values in case there are no tracks in the event
    m_ft1energy = CTBBestEnergy;
    if( m_ft1energy==0) m_ft1energy = EvtEnergyCorr;

    m_ft1theta = 666; m_ft1phi = 666; m_ft1ra   = 666;
    m_ft1dec   = 666; m_ft1zen = 666; m_ft1azim = 666;
    m_ft1l = 666; m_ft1b = 666;
    m_ft1convpointx = 999; m_ft1convpointy = 999; m_ft1convpointz = 999; 
    m_ft1convlayer = -1;
    m_ft1livetime = -1;
    m_ft1livetime = EvtLiveTime;

    if( TkrNumTracks==0) return;

    CLHEP::Hep3Vector convPoint;
    CLHEP::Hep3Vector glastDir;
    if( CTBBestZDir!=0){ // check that this was set
        glastDir= CLHEP::Hep3Vector(CTBBestXDir, CTBBestYDir, CTBBestZDir);
    }else{
        glastDir= CLHEP::Hep3Vector(Tkr1XDir, Tkr1YDir, Tkr1ZDir);
    }

    m_ft1convlayer   = Tkr1FirstLayer;

    glastDir = - glastDir.unit();

    if( exp==0 ) return; // do not fill in orientation-dependent stuff if no pointing info

    // celestial coords in degree
    std::map<std::string,double> cel_coords = getCelestialCoords(*exp, glastDir);
    m_ft1ra   = cel_coords["RA"];
    m_ft1dec  = cel_coords["DEC"];
    astro::SkyDir my_dir(m_ft1ra, m_ft1dec);
    m_ft1l = my_dir.l();
    m_ft1b = my_dir.b();
    m_ft1zen  = cel_coords["ZENITH_THETA"];
    m_ft1azim = cel_coords["EARTH_AZIMUTH"];

    // instrument coords in degree
    m_ft1theta = glastDir.theta()*180/M_PI;
    double phi_deg = glastDir.phi(); 
    if( phi_deg<0 ) phi_deg += 2*M_PI;
    m_ft1phi =  phi_deg*180/M_PI;

}
//------------------------------------------------------------------------------
std::map<std::string, double> 
FT1worker::getCelestialCoords(const Event::Exposure& exp, const CLHEP::Hep3Vector glastDir)
{
    using namespace astro;
    using CLHEP::Hep3Vector;

    std::map<std::string, double> fields;


    // create a transformation object -- first get local directions
    SkyDir zsky( exp.RAZ(), exp.DECZ() );
    SkyDir xsky( exp.RAX(), exp.DECX() );
    // orthogonalize, since interpolation and transformations destory orthogonality (limit is 10E-8)
    Hep3Vector xhat = xsky() -  xsky().dot(zsky()) * zsky() ;
    PointingTransform toSky( zsky, xhat );

    // make zenith (except for oblateness correction) unit vector
    Hep3Vector position( exp.posX(),  exp.posY(),  exp.posZ() );
    SkyDir zenith(position.unit());

    SkyDir sdir = toSky.gDir(glastDir);

    //zenith_theta and earth_azimuth
    double zenith_theta = sdir.difference(zenith); 
    if( fabs(zenith_theta)<1e-8) zenith_theta=0;

    // all this to do the azimuth angle :-(
    Hep3Vector north_pole(0,0,1);
    Hep3Vector east_dir( north_pole.cross(zenith()).unit() ); // east is perp to north_pole and zenith
    Hep3Vector north_dir( zenith().cross(east_dir));

    double earth_azimuth=atan2( sdir().dot(east_dir), sdir().dot(north_dir) );
    if( earth_azimuth <0) earth_azimuth += 2*M_PI; // to 0-360 deg.
    if( fabs(earth_azimuth)<1e-8) earth_azimuth=0;

    fields["RA"]            = sdir.ra();
    fields["DEC"]           = sdir.dec();
    fields["ZENITH_THETA"]  = zenith_theta*180/M_PI;
    fields["EARTH_AZIMUTH"] = earth_azimuth*180/M_PI;

    return fields;
}


