/** @file PtValsAlg.cxx
@brief declaration and definition of the class PtValsAlg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/PtValsAlg.cxx,v 1.11 2008/10/01 07:02:16 lsrea Exp $

*/

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/SmartDataPtr.h"

// Event for access to time
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

// stuff for exposure info (not really MC)
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/Exposure.h"

// to write a Tree with pointing info
#include "ntupleWriterSvc/INTupleWriterSvc.h"


#include "facilities/Util.h"
#include "facilities/Timestamp.h"

#include "astro/JulianDate.h"
#include "astro/GPS.h"

//
#include "AnalysisNtuple/PointingInfo.h"

namespace { // anonymous namespace for file-global
    astro::GPS* gps(0);  // pointer to relevant GPS entry
}

using namespace AnalysisNtuple;

/** 
* \class PtValsAlg
*
* \brief This is an Algorithm designed to store pointing information in the tuple
* \author Toby Burnett
* 
*/


class PtValsAlg : public Algorithm {
public:
    PtValsAlg(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();


private: 
    PointingInfo m_pointingInfo;

    StringProperty m_root_tree;
    StringArrayProperty m_pointingHistory;///< history file name and launch date

    INTupleWriterSvc* m_rootTupleSvc;
    IDataProviderSvc* m_pEventSvc;

    astro::PointingHistory* m_history;
    bool m_horizontal;
    bool m_fillNtuple;

};
//------------------------------------------------------------------------


static const AlgFactory<PtValsAlg>  Factory;
const IAlgFactory& PtValsAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
PtValsAlg::PtValsAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator) , m_horizontal(false)
{
    // declare properties with setProperties calls

    declareProperty("pointing_info_tree_name",  m_root_tree="MeritTuple");
    // doublet, filename and launch date
    declareProperty("PointingHistory",   m_pointingHistory); 
    declareProperty("FillNtuple",        m_fillNtuple=true);

}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode PtValsAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get a pointer to RootTupleSvc
    if( (service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " RootTupleSvc is not available" << endreq;
        m_rootTupleSvc=0;
        sc = StatusCode::FAILURE;
    } else if( !m_root_tree.value().empty() ) {   
     
          if(m_fillNtuple) m_pointingInfo.setPtTuple(m_rootTupleSvc, m_root_tree.value()); 
    }


    // get the GPS instance: 

    gps = astro::GPS::instance();

    //set the input file to be used as the pointing database, if used
    if( m_pointingHistory.value().empty()){
        log << MSG::WARNING << "No history file specified, using default" << endreq;
    }else{
        std::string filename(m_pointingHistory.value()[0]);
        facilities::Util::expandEnvVar(&filename);
        double offset = 0;
        std::string jStr;
        if( m_pointingHistory.value().size()>1){
            std::string field(m_pointingHistory.value()[1]);
            if(! field.empty() ) { // allow null string
                facilities::Timestamp jt(m_pointingHistory.value()[1]);
                offset = (astro::JulianDate(jt.getJulian())
                    - astro::JulianDate::missionStart())
                    *astro::JulianDate::secondsPerDay;
                astro::JulianDate jDate(astro::JulianDate(jt.getJulian()));
                jStr = jDate.getGregorianDate();
            }
        }

        if( m_pointingHistory.value().size()>2){
            std::string field(m_pointingHistory.value()[2]);
            m_horizontal =! field.empty();
        }
        log << MSG::INFO << "Loading Pointing History File : " << filename <<endreq;
        if( offset>0 ){
            log << MSG::INFO   << " with MET offset " ;
            log.precision(12);
            log << offset << " ";
            log.precision(6);
            log << jStr << endreq;
        }
        if( m_horizontal){
            log << MSG::INFO 
                << "   Will override x-direction to be horizontal"<<endreq;
        }
        gps->setPointingHistoryFile(filename, offset, m_horizontal);
    }



    IDataProviderSvc* eventsvc = 0;
    sc = serviceLocator()->service( "EventDataSvc", eventsvc, true );
    if(sc.isFailure()){
        log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
        return sc;
    }
    m_pEventSvc = eventsvc;
   
    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode PtValsAlg::execute()
{
    //StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    //
    // Purpose: set tuple items
 
   SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);

   // get event time from header or merit 
   // and look up position info from the history  

   double etime;

   if(header==0) {
       void* ptr;
       m_rootTupleSvc->getItem(m_root_tree,"PtTime",ptr);
       etime = *reinterpret_cast<double*>(ptr);
   } else {
       etime = header->time();
   }
    gps->time(etime);

    // and create the tuple
    if (m_fillNtuple) m_pointingInfo.execute( *gps );
    
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode PtValsAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    return sc;
}



