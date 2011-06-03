/** @file FT1Alg.cxx
@brief Declaration and implementation of Gaudi algorithm FT1Alg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/FT1Alg.cxx,v 1.29 2011/05/31 12:19:12 bregeon Exp $
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

    const double R2D = 180./M_PI;
    const int _invalidEventClass = -101;

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
    Item EvtElapsedTime;
    Item EvtLiveTime;
    Item EvtEnergyCorr;
    Item TkrNumTracks;
    Item VtxXDir, VtxYDir, VtxZDir;
    Item VtxX0, VtxY0, VtxZ0;
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
    Item CTBParticleType;
    Item Cal1MomXDir, Cal1MomYDir, Cal1MomZDir;

    //FT1 entries to create
    unsigned int m_ft1eventid;
    float m_ft1energy;
    float m_ft1theta,m_ft1phi,m_ft1ra,m_ft1dec,m_ft1l,m_ft1b;
    float m_ft1zen,m_ft1azim;
    double m_ft1livetime;
    //float m_ft1convpointx,m_ft1convpointy,m_ft1convpointz,
    float m_ft1convlayer;
    int   m_ft1eventclass;
    float m_ft1calzen,m_ft1calazim,m_ft1calra,m_ft1caldec,m_ft1call,m_ft1calb;
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
            // explanation from Joanne for the following: 
            //"The type of (*r)  is a  3-element array of double, so applying a subscript
            // should give you one of the doubles in that array. r itself is a pointer to something 3 doubles long "
            double x((*r)[0]), y((*r)[1]), z((*r)[2]);
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
, CTBParticleType("CTBParticleType")
, Cal1MomXDir("Cal1MomXDir")
, Cal1MomYDir("Cal1MomYDir")
, Cal1MomZDir("Cal1MomZDir")

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
    //addItem( "FT1ConvPointX",    m_ft1convpointx); //gone
    //addItem( "FT1ConvPointY",    m_ft1convpointy);
    //addItem( "FT1ConvPointZ",    m_ft1convpointz);
    addItem( "FT1ConvLayer",     m_ft1convlayer);
    addItem( "FT1Livetime",      m_ft1livetime);
    addItem( "FT1EventClass",    m_ft1eventclass);
    addItem( "FT1CalZenithTheta",  m_ft1calzen);
    addItem( "FT1CalEarthAzimuth", m_ft1calazim);
    addItem( "FT1CalRa",           m_ft1calra);
    addItem( "FT1CalDec",          m_ft1caldec);
    addItem( "FT1CalL",            m_ft1call);
    addItem( "FT1CalB",            m_ft1calb);
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
<td>D<td>  (s) Cumulative live time, from start of run, or mission <b>do not use</b>
<tr><td> FT1ConvLayer 
<td>F<td>  Starting layer of the best track found in the tracker 
(Layer 0 is the one closest to the calorimeter.)
<tr><td> FT1ConvPoint[X/Y/Z] 
<td>F<td>  REMOVED! <b>Do not use; no longer filled!</b>
<tr><td> FT1EventClass
<td>I<td> 
   Combination of CTBParticleType and CTBClassLevel
   see: https://confluence.slac.stanford.edu/display/SCIGRPS/CTBClassLevel+Definition+for+P7
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
    //m_ft1convpointx = 999; m_ft1convpointy = 999; m_ft1convpointz = 999; //gone 
    m_ft1convlayer = -1;
    m_ft1livetime = -1;
    m_ft1livetime = EvtLiveTime;
    m_ft1eventclass = _invalidEventClass;
    m_ft1calra  = 666; m_ft1caldec  = 666;
    m_ft1calzen = 666; m_ft1calazim = 666;
    m_ft1call   = 666; m_ft1calb    = 666;

    // first calculate the EventClass, pass-7 style

    int particleType = (int)floor(CTBParticleType + 0.001);
    int classLevel   = (int)floor(100*CTBClassLevel + 0.001);

    if (particleType==0) { // some kind of gamma

    // The code below replaces this:
        // if     (classLevel == 5)   {m_ft1eventclass = 1;}   // 0.05
        // else if(classLevel == 10)  {m_ft1eventclass = 2;}   // 0.1
        // else if(classLevel == 20)  {m_ft1eventclass = 3;}   // 0.2
        // else if(classLevel == 30)  {m_ft1eventclass = 4;}   // 0.3
        // else if(classLevel == 50)  {m_ft1eventclass = 5;}   // 0.5
        // else if(classLevel == 100) {m_ft1eventclass = 6;}   // 1.0
        // else if(classLevel == 200) {m_ft1eventclass = 7;}   // 2.0
        // else if(classLevel == 300) {m_ft1eventclass = 8;}   // 3.0
        // else if(classLevel == 400) {m_ft1eventclass = 9;}   // 4.0
        // else if(classLevel == 500) {m_ft1eventclass = 10;}  // 5.0        
 
        switch (classLevel) {
            // these are the only allowed values
            case   5: m_ft1eventclass =  1; break;  // 0.05
            case  10: m_ft1eventclass =  2; break;  // 0.1
            case  20: m_ft1eventclass =  3; break;  // 0.2
            case  30: m_ft1eventclass =  4; break;  // 0.3
            case  50: m_ft1eventclass =  5; break;  // 0.5
            case 100: m_ft1eventclass =  6; break;  // 1.0
            case 200: m_ft1eventclass =  7; break;  // 2.0
            case 300: m_ft1eventclass =  8; break;  // 3.0
            case 400: m_ft1eventclass =  9; break;  // 4.0
            case 500: m_ft1eventclass = 10; break; // 5.0
            default:  m_ft1eventclass = _invalidEventClass;  // should never happen   
        }            
    } else if (particleType>=-1 && particleType<=4){
        if(classLevel==100||classLevel==200||classLevel==300 ||
            (particleType==-1&&classLevel==400)) { // only levels for these types
            m_ft1eventclass = particleType*100 + classLevel/100;
        }
    }

    //
    // Projection of tracker and calorimeter direction into the sky and the Earth
    //
    
    // First deal with the tracker projection
    // although CTBBest[X/Y/Z]Dir can be set even if there's no track, 
    // by Cal-only events, for example!
    Hep3Vector glastDir;
    // if there is a BestDir use it
    if(CTBBestZDir!=0)
        glastDir= Hep3Vector(CTBBestXDir, CTBBestYDir, CTBBestZDir);
    // if there is no BestDir but a Track exists, use it
    else if(TkrNumTracks>0)
        glastDir= Hep3Vector(Tkr1XDir, Tkr1YDir, Tkr1ZDir);
        
    // if we have a tracker direction to project, do it now
    if(glastDir.mag2() != 0) {
        // instrument coords
        m_ft1convlayer   = Tkr1FirstLayer;

        m_ft1theta = (-glastDir).theta()*R2D;
        double phi_deg = (-glastDir).phi(); 
        if( phi_deg<0 ) phi_deg += 2*M_PI;
        m_ft1phi =  phi_deg*R2D;

        // celestial coordinates
        // glastDir points down... 
        // toSky converts a *particle* direction
        // into a direction on the sky, so the minus-sign is taken care of
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
        // east is perp to north_pole and zenith
        Hep3Vector east_dir( north_pole.cross(zenith()).unit() );
        Hep3Vector north_dir( zenith().cross(east_dir));

        double earth_azimuth=atan2( sdir().dot(east_dir), sdir().dot(north_dir) );
        if( earth_azimuth <0) earth_azimuth += 2*M_PI; // to 0-360 deg.
        if( fabs(earth_azimuth)<1e-8) earth_azimuth=0;
        m_ft1zen  = zenith_theta*R2D;;
        m_ft1azim = earth_azimuth*R2D;
    } // done with tracker projections
    
    
    // Now deal with the calorimater projection
    // for now we project only if we have results from the moment analysis
    Hep3Vector glastCalDir;
    if ( Cal1MomZDir != 0 ) {
        glastCalDir = Hep3Vector(Cal1MomXDir, Cal1MomYDir, Cal1MomZDir);
      
        // celestial coordinates
        // The Cal axis points up, so we *do* need the minus sign, here.
        SkyDir scaldir( gps->toSky(-glastCalDir) ); 
        m_ft1calra  = scaldir.ra();
        m_ft1caldec = scaldir.dec();
        m_ft1call   = scaldir.l();
        m_ft1calb   = scaldir.b();

        // local Zenith coordinates
        SkyDir zenith(gps->zenithDir());  // pointing direction
        double calzenith_theta = scaldir.difference(zenith); 
        if( fabs(calzenith_theta)<1e-8) calzenith_theta=0;

        // all this to do the azimuth angle :-(
        Hep3Vector north_pole(0,0,1);
        // east is perp to north_pole and zenith
        Hep3Vector east_dir( north_pole.cross(zenith()).unit() );
        Hep3Vector north_dir( zenith().cross(east_dir));

        double calearth_azimuth=atan2( scaldir().dot(east_dir), scaldir().dot(north_dir) );
        if( calearth_azimuth <0) calearth_azimuth += 2*M_PI; // to 0-360 deg.
        if( fabs(calearth_azimuth)<1e-8) calearth_azimuth=0;

        m_ft1calzen  = calzenith_theta*R2D;
        m_ft1calazim = calearth_azimuth*R2D;      
   } // done with calorimeter projections

 }
