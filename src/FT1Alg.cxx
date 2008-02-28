/** @file FT1Alg.cxx
@brief Declaration and implementation of Gaudi algorithm FT1Alg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/FT1Alg.cxx,v 1.7 2008/02/27 20:11:50 burnett Exp $
*/
// Include files

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"


// Event for access to time
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "astro/GPS.h"

#include "FluxSvc/IFluxSvc.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <cassert>
#include <map>

// forward declatation of the worker
class FT1worker;

namespace { // anonymous namespace for file-global
    unsigned int nbOfEvtsInFile(100000);
    std::string treename("MeritTuple");

    astro::GPS* gps(0);  // pointer to relevant GPS entry
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
    IDataProviderSvc* m_pEventSvc;

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static const AlgFactory<FT1Alg>  Factory;
const IAlgFactory& FT1AlgFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class FT1worker{ 
public:
    FT1worker();

    void evaluate();
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
FT1Alg::FT1Alg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
, m_count(0)
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

    IDataProviderSvc* eventsvc = 0;
    sc = serviceLocator()->service( "EventDataSvc", eventsvc, true );
    if(sc.isFailure()){
        log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
        return sc;
    }
    m_pEventSvc = eventsvc;


    // get a pointer to RootTupleSvc 
    if( (sc = service("RootTupleSvc", rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
        return sc;
    }

    m_worker =  new FT1worker();

    // get the GPS instance: either from FluxSvc or local, non-MC mode
    IFluxSvc* fluxSvc(0);
    if( service("FluxSvc", fluxSvc, true).isFailure() ){
        gps = astro::GPS::instance();
    }else{
        gps = fluxSvc->GPSinstance();
    }
    return sc;
}

StatusCode FT1Alg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    m_count++;
   SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);

    // get event time from header, and make sure gps is in sync 
    double etime(header->time());
    gps->time(etime); 

    m_worker->evaluate();
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
#if 0
void FT1worker::evaluate(const Event::Exposure* exp)
#else
void FT1worker::evaluate()
#endif
{
    using CLHEP::Hep3Vector;
    using astro::SkyDir;

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
    m_ft1convlayer   = Tkr1FirstLayer;

    Hep3Vector glastDir;
    if( CTBBestZDir!=0){ // check that this was set
        glastDir= Hep3Vector(CTBBestXDir, CTBBestYDir, CTBBestZDir);
    }else{
        glastDir= Hep3Vector(Tkr1XDir, Tkr1YDir, Tkr1ZDir);
    }

    // instrument coords

    m_ft1theta = (-glastDir).theta()*180/M_PI;
    double phi_deg = (-glastDir).phi(); 
    if( phi_deg<0 ) phi_deg += 2*M_PI;
    m_ft1phi =  phi_deg*180/M_PI;

    // celestial coordinates

    // transform 
    SkyDir sdir( gps->toSky(glastDir) ); 
    m_ft1ra  = sdir.ra();
    m_ft1dec = sdir.dec();
    m_ft1l   = sdir.l();
    m_ft1b   = sdir.b();

    // local Zenith coordinates

    SkyDir zenith(gps->zenithDir());  // pointing direction
    double zenith_theta = sdir.difference(zenith); 
    if( fabs(zenith_theta)<1e-8) zenith_theta=0;

    // all this to do the azimuth angle :-(
    Hep3Vector north_pole(0,0,1);
    Hep3Vector east_dir( north_pole.cross(zenith()).unit() ); // east is perp to north_pole and zenith
    Hep3Vector north_dir( zenith().cross(east_dir));

    double earth_azimuth=atan2( sdir().dot(east_dir), sdir().dot(north_dir) );
    if( earth_azimuth <0) earth_azimuth += 2*M_PI; // to 0-360 deg.
    if( fabs(earth_azimuth)<1e-8) earth_azimuth=0;

    m_ft1zen  = zenith_theta*180/M_PI;;
    m_ft1azim = earth_azimuth*180/M_PI;

}

