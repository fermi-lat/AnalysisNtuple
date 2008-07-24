/** @file FT1Alg.cxx
@brief Declaration and implementation of Gaudi algorithm FT1Alg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/FT1Alg.cxx,v 1.16 2008/07/24 03:52:54 burnett Exp $
*/
// Include files

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "CalibData/Nas/CalibLATAlignment.h"

// Event for access to time
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "astro/GPS.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include "CalibData/Nas/CalibLATAlignment.h" // defines interface for calibration object 
#include "CalibSvc/ICalibPathSvc.h" // helpful service for forming calib TDS paths

//#include "CalibSvc/ICalibDataSvc.h"

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

    IDataProviderSvc* m_pCalibDataSvc;
    ICalibPathSvc*    m_pCalibPathSvc;

    BooleanProperty m_aberrate;  ///< set true to enable aberration correction
    StringProperty m_flavor;     ///< set to a flavor to enable the corrrection

    std::string m_path;    ///< calibration service
    

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
// LSR 14-Jul-08 code for ntuple types; potential changes here!
    TypedItem<unsigned int, 'i'> EvtRun; 
    TypedItem<unsigned long long , 'l'> EvtEventId64;

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
    Item CTBBestEnergyProb;
    Item CTBBestEnergyRatio;
    Item CTBCORE;
    Item CTBClassLevel;
    

    //FT1 entries to create
    unsigned int m_ft1eventid;
    float m_ft1energy;
    float m_ft1theta,m_ft1phi,m_ft1ra,m_ft1dec,m_ft1l,m_ft1b;
    float m_ft1zen,m_ft1azim;
    double m_ft1livetime;
    float m_ft1convpointx,m_ft1convpointy,m_ft1convpointz,m_ft1convlayer;
    float m_ft1eventclass;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FT1Alg::FT1Alg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
, m_count(0)
{
    declareProperty("TreeName",  treename="MeritTuple");
    declareProperty("NbOfEvtsInFile", nbOfEvtsInFile=100000);
    declareProperty("CorrectForAberration", m_aberrate=false);
    declareProperty("AlignmentFlavor"     , m_flavor="");
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

    // get the GPS instance: 
    gps = astro::GPS::instance();

    // enable aberration correction if requested
    if( m_aberrate.value() ){
        gps->enableAberration(true);
        log << MSG::INFO << "Correction for aberration is enabled" << endreq;
    }
    // stuff for getting the calibration -- this code grateful curtesy of J. Bogart!
    if( !m_flavor.value().empty() ) {
        if( (sc=service("CalibDataSvc", m_pCalibDataSvc, true)).isFailure()) { // needed for any access to calib. tds
            log << MSG::ERROR << " failed to get the CalibDataSvc" << endreq;
            return sc;
        }
        if( (sc=service("CalibDataSvc", m_pCalibPathSvc, true)).isFailure()){ // preferred way to form paths into calib tds
            log << MSG::ERROR << " failed to get the CalibDataSvc (or maybe CalibPathSvc)" << endreq;
            return sc;
        }

        m_path = m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_NAS_LATAlignment, m_flavor); // form path
        log << MSG::INFO << "Using alignment flavor "<< m_flavor.value() << endreq;

    }
    return sc;
}


StatusCode FT1Alg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    m_count++;

    if(! m_flavor.value().empty() ){
        // update alignment if necessary -- this code grateful curtesy of J. Bogart!


        static unsigned serial = 0; // serial number of last calibration used

        SmartDataPtr<CalibData::CalibLATAlignment> alignCalib(m_pCalibDataSvc, m_path);
        if ( !alignCalib ) {
            log << MSG::ERROR << "Failed access to CalibLATAlignment via smart ptr with path"
                << m_path << endreq;
            return StatusCode::FAILURE;
        }
        unsigned newSerial = alignCalib->getSerNo();

        if( serial != newSerial) {
            serial = newSerial;

            static double arcsec2deg(M_PI/648000);
            const CalibData::ALIGN_ROT* r = alignCalib->getR();
            double x(*r[0]), y(*r[1]), z(*r[2]);
            gps->setAlignmentRotation( 
                CLHEP::HepRotationX(x * arcsec2deg)
                *CLHEP::HepRotationY(y* arcsec2deg)
                *CLHEP::HepRotationZ(z* arcsec2deg)
                );
            log << MSG::INFO << "setting alignment parameters: ("
                << x<<", "<< y<< "," << z << ") arcsec" << endreq;
        }
    }

    
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
, EvtEventId64("EvtEventId64")
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
, CTBBestEnergyProb("CTBBestEnergyProb")
, CTBBestEnergyRatio("CTBBestEnergyRatio")
, CTBCORE("CTBCORE")
, CTBClassLevel("CTBClassLevel")

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
    addItem( "FT1EventClass",    m_ft1eventclass);
}

/** @page anatup_vars 

@section FT1  FT1 Variables

see <a href="http://glast.gsfc.nasa.gov/ssc/dev/fits_def/definitionFT1.html">FT1 definition</a>

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> FT1EventId  
<td>U<td>  RunNumber*100000 + EventNumber (overflows for real data!)
<tr><td> FT1Energy   
<td>F<td>  (MeV) estimate for energy  
<tr><td> FT1Theta, FT1Phi 
<td>F<td>  (deg) reconstructed direction with respect to instrument coordinate system      
<tr><td> FT1Ra, FT1Dec 
<td>F<td>  (deg) reconstructed direction in equatorial coordinates       
<tr><td> FT1ZenithTheta
<td>F<td>  (deg) reconstucted angle with respect to local zenith
<tr><td> FT1EarthAzimuth 
<td>F<td>  (deg) reconstucted heading: 0 is North, 90 East, etc.
<tr><td> FT1L, FT1B 
<td>F<td>  (deg) galactic longitude and latitude of reconstructed direction
<tr><td> FT1Livetime 
<td>D<td>  (s) Cumulative live time, from start of run, or mission <b>do not use<b>
<tr><td> FT1ConvLayer 
<td>F<td>  Starting layer of the best track found in the tracker 
(Layer 0 is the one closest to the calorimeter.)
<tr><td> FT1ConvPoint[X/Y/Z] 
<td>F<td>  <b>Do not use; no longer filled!</b>
<tr><td> FT1EventClass
<td>F<td> 
   CTBBestEnergy<=10 || CTBBestEnergyRatio>=5 -> 0 <br>
   else CTBClassLevel>2 && CTBCORE>0. && CTBBestEnergyProb>0.1 -> 3<br>
   else CTBClassLevel>1 && CTBCORE>0.1 && CTBBestEnergyProb>0.1 -> 2<br>
   else CTBClassLevel>0 && CTBCORE>0. && CTBBestEnergyProb>0. -> 1<br>
   else 0
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
    m_ft1eventid =  nbOfEvtsInFile *EvtRun.value()  + EvtEventId64.value();

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
    // glastDir points down... 
    // toSky converts a *particle* direction
    // into a direction on the sky, so the minus-sign is taken care of!
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

    // more useful class level

    m_ft1eventclass = 0;

    if ( (CTBBestEnergy<=10) || (CTBBestEnergyRatio>=5) ) { // Apply minimal cut first.
        m_ft1eventclass = 0;
    } else if ( (CTBClassLevel>2) && (CTBCORE>0.1) && (CTBBestEnergyProb>0.1) ) {
        m_ft1eventclass = 3;
    } else if ( (CTBClassLevel>1) && (CTBCORE>0.1) && (CTBBestEnergyProb>0.1) ) {
        m_ft1eventclass = 2;
    } else if ( (CTBClassLevel>0) && (CTBCORE>0.) && (CTBBestEnergyProb>0.) ) {
        m_ft1eventclass = 1;
    }

}
