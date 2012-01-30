/** @file FT1Alg.cxx
@brief Declaration and implementation of Gaudi algorithm FT1Alg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/FT1Alg.cxx,v 1.27.10.5 2011/10/27 02:17:16 heather Exp $
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

#include "facilities/Util.h"

#include "evtUtils/EventClass.h"


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
 
    /// this guy handles the event class stuff
    evtUtils::EventClass* m_evtClass;


    //counter
    int m_count;
    IDataProviderSvc* m_pEventSvc;

    IDataProviderSvc* m_pCalibDataSvc;
    ICalibPathSvc*    m_pCalibPathSvc;

    INTupleWriterSvc* m_rootTupleSvc;
   

    BooleanProperty m_aberrate;  ///< set true to enable aberration correction
    StringProperty m_flavor;     ///< set to a flavor to enable the corrrection
    StringProperty m_classDefs;  ///< xml file with events class defs

    std::string m_path;    ///< calibration service

    std::string m_defaultEvtClassDefFile; ///< default value for input XML
    

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//static const AlgFactory<FT1Alg>  Factory;
//const IAlgFactory& FT1AlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(FT1Alg);

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
    Item Tkr1XDir, Tkr1YDir, Tkr1ZDir;
    Item Tkr1FirstLayer;
    Item CTBBestEnergy;
    Item CTBBestXDir;
    Item CTBBestYDir;
    Item CTBBestZDir;

    //FT1 entries to create
    unsigned int m_ft1eventid;
    float m_ft1energy;
    float m_ft1theta,m_ft1phi,m_ft1ra,m_ft1dec,m_ft1l,m_ft1b;
    float m_ft1zen,m_ft1azim;
    double m_ft1livetime;
    float m_ft1convlayer;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FT1Alg::FT1Alg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
, m_count(0), m_defaultEvtClassDefFile("$(EVTUTILSXMLPATH)/EvtClassDefs_P7V6.xml")
{
    declareProperty("TreeName",  treename="MeritTuple");
    declareProperty("NbOfEvtsInFile", nbOfEvtsInFile=100000);
    declareProperty("CorrectForAberration", m_aberrate=false);
    declareProperty("AlignmentFlavor"     , m_flavor="");
    facilities::Util::expandEnvVar(&m_defaultEvtClassDefFile);
    declareProperty("EventClassDefinitions", m_classDefs=m_defaultEvtClassDefFile);
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
    if( (sc = service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
        return sc;
    }
    rootTupleSvc = m_rootTupleSvc;

    m_worker =  new FT1worker;

    std::string inputEvtClassDefFile = m_classDefs;
    facilities::Util::expandEnvVar(&inputEvtClassDefFile);
    if (inputEvtClassDefFile.compare(m_defaultEvtClassDefFile) != 0)
        log << MSG::WARNING << "Event Class Definition File " <<
            inputEvtClassDefFile << " does not match the default " 
            << m_defaultEvtClassDefFile << endreq;
    
    m_evtClass = evtUtils::EventClass::loadFromXml(m_classDefs);
    if ( m_evtClass == 0 ) {
        log << MSG::ERROR << " failed load event class definitions from " << m_classDefs << endreq;
        return StatusCode::FAILURE;      
    }
    std::list<std::string> ft1EvtClassNames;
    m_evtClass->getEvtMapNames(ft1EvtClassNames);

    for ( std::list<std::string>::const_iterator itr = ft1EvtClassNames.begin(); itr != ft1EvtClassNames.end(); itr++) {
        addItem(itr->c_str(),*(m_evtClass->getShortMapPtr(*itr)));
    }

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

    // We need the tree and we need to jump through some hoops to get it
    void* treePtrVoid(0);
    m_rootTupleSvc->getItem("MeritTuple","",treePtrVoid);
    TTree* treePtr = static_cast<TTree*>(treePtrVoid);
    if ( treePtr == 0 ) {
        log << MSG::ERROR << " Failed to get MeritTuple to build class defintions." << m_classDefs << endreq;
        return StatusCode::FAILURE;            
    }

    // This next step doesn't do anything if the treePtr matches the cached value
    if ( ! m_evtClass->initializeShortCuts(*treePtr) ) {
        log << MSG::ERROR << " Failed to initialize Event Class definitions against merit tuple." << m_classDefs << endreq;
        return StatusCode::FAILURE;                 
    }

    // now set the bits
    if ( ! m_evtClass->fillShortCutMaps() ) {
        log << MSG::ERROR << " Failed to fill cut maps." << endreq;
        return StatusCode::FAILURE;                      
    }

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
, Tkr1XDir("Tkr1XDir")
, Tkr1YDir("Tkr1YDir")
, Tkr1ZDir("Tkr1ZDir")
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
    addItem( "FT1ConvLayer",     m_ft1convlayer);
    addItem( "FT1Livetime",      m_ft1livetime);
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


    // CTBBest[X/Y/Z]Dir can be set even if there's no track, 
    // by Cal-only events, for example!

    Hep3Vector glastDir;
    if( CTBBestZDir!=0){ // check that this was set
        glastDir= Hep3Vector(CTBBestXDir, CTBBestYDir, CTBBestZDir);
    } else {
        if(TkrNumTracks>0) {
            // Track exists, use it
            glastDir= Hep3Vector(Tkr1XDir, Tkr1YDir, Tkr1ZDir);
        } else {
            // nothing to do, leave quietly!
            return;
        }
    }

    // instrument coords

    m_ft1convlayer   = Tkr1FirstLayer;

    m_ft1theta = (-glastDir).theta()*R2D;
    double phi_deg = (-glastDir).phi(); 
    if( phi_deg<0 ) phi_deg += 2*M_PI;
    m_ft1phi =  phi_deg*R2D;

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

    m_ft1zen  = zenith_theta*R2D;;
    m_ft1azim = earth_azimuth*R2D;

 }
