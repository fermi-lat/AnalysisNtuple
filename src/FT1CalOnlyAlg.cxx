/** @file FT1CalOnlyAlg.cxx
@brief Declaration and implementation of Gaudi algorithm FT1CalOnlyAlg

$Header:  $
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

#include <cassert>
#include <map>

// forward declaration of the worker
class FT1CalOnlyWorker;

namespace { // anonymous namespace for file-global
    unsigned int nbOfEvtsInFile(100000);
    std::string treename("MeritTuple");

    const double R2D = 180./M_PI;
    const int _invalidEventClass = -101;

    astro::GPS* gps(0);  // pointer to relevant GPS entry

#include "Item.h"
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class FT1CalOnlyAlg
@brief Stripped-down version of the TF1Alg class to be run when only the CAL variables
are defined.

If you wonder why one would be interested in this, well... it was originally written
for skimming in two steps the GC events based on the calorimeter direction.
*/
class FT1CalOnlyAlg : public Algorithm {

public:
    FT1CalOnlyAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
private:
    /// this guy does the work!
    FT1CalOnlyWorker * m_worker;
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

DECLARE_ALGORITHM_FACTORY(FT1CalOnlyAlg);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class FT1CalOnlyWorker{ 
public:
  FT1CalOnlyWorker();
  
  void evaluate();
  // these all float or double
  Item Cal1MomXDir, Cal1MomYDir, Cal1MomZDir;

  //FT1 entries to create
  float m_ft1calzen,m_ft1calazim,m_ft1calra,m_ft1caldec,m_ft1call,m_ft1calb;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FT1CalOnlyAlg::FT1CalOnlyAlg(const std::string& name, ISvcLocator* pSvcLocator) 
  : Algorithm(name, pSvcLocator)
  , m_count(0)
{
  declareProperty("TreeName", treename="MeritTuple");
  declareProperty("NbOfEvtsInFile", nbOfEvtsInFile=100000);
  declareProperty("CorrectForAberration", m_aberrate=false);
  declareProperty("AlignmentFlavor", m_flavor="");
}

StatusCode FT1CalOnlyAlg::initialize()
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
  
  m_worker =  new FT1CalOnlyWorker();
  
  // get the GPS instance: 
  gps = astro::GPS::instance();
  
  // enable aberration correction if requested
  if( m_aberrate.value() ){
    gps->enableAberration(true);
    log << MSG::INFO << "Correction for aberration is enabled" << endreq;
  }
  // stuff for getting the calibration -- this code grateful curtesy of J. Bogart!
  if( !m_flavor.value().empty() ) {
    if( (sc=service("CalibDataSvc", m_pCalibDataSvc, true)).isFailure()) {
      // needed for any access to calib. tds
      log << MSG::ERROR << " failed to get the CalibDataSvc" << endreq;
      return sc;
    }
    if( (sc=service("CalibDataSvc", m_pCalibPathSvc, true)).isFailure()){
      // preferred way to form paths into calib tds
      log << MSG::ERROR << " failed to get the CalibDataSvc (or maybe CalibPathSvc)" <<
        endreq;
      return sc;
    }
    
    m_path = m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_NAS_LATAlignment, m_flavor);
    log << MSG::INFO << "Using alignment flavor "<< m_flavor.value() << endreq;
    
  }
  return sc;
}


StatusCode FT1CalOnlyAlg::execute()
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
      // "The type of (*r)  is a  3-element array of double, so applying a subscript
      // should give you one of the doubles in that array. r itself is a pointer to
      // something 3 doubles long "
      double x((*r)[0]), y((*r)[1]), z((*r)[2]);
      gps->setAlignmentRotation(CLHEP::HepRotationX(x * arcsec2deg)
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

StatusCode FT1CalOnlyAlg::finalize()
{
  return StatusCode::SUCCESS;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FT1CalOnlyWorker::FT1CalOnlyWorker()
  // initialize pointers to current items
  :Cal1MomXDir("Cal1MomXDir")
  ,Cal1MomYDir("Cal1MomYDir")
  ,Cal1MomZDir("Cal1MomZDir")
{
  //now create new items 
  addItem( "FT1CalZenithTheta",  m_ft1calzen);
  addItem( "FT1CalEarthAzimuth", m_ft1calazim);
  addItem( "FT1CalRa",           m_ft1calra);
  addItem( "FT1CalDec",          m_ft1caldec);
  addItem( "FT1CalL",            m_ft1call);
  addItem( "FT1CalB",            m_ft1calb);
}


#if 0
void FT1CalOnlyWorker::evaluate(const Event::Exposure* exp)
#else
  void FT1CalOnlyWorker::evaluate()
#endif
{
  using CLHEP::Hep3Vector;
  using astro::SkyDir;

  m_ft1calra  = 666; m_ft1caldec  = 666;
  m_ft1calzen = 666; m_ft1calazim = 666;
  m_ft1call   = 666; m_ft1calb    = 666;
  
  //
  // Projection of calorimeter direction into the sky and the Earth
  //
  

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
  }
  
}
