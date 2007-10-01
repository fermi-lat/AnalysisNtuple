/** @file PtValsAlg.cxx
@brief declaration and definition of the class PtValsAlg

$Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/PtValsAlg.cxx,v 1.5 2006/11/11 22:01:59 burnett Exp $

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
#include "astro/PointingHistory.h"

//
#include "AnalysisNtuple/PointingInfo.h"

using namespace AnalysisNtuple;


/** 
* \class PtValsAlg
*
* \brief This is an Algorithm designed to store pointing information in the tuple
* \author Toby Burnett
* 
* $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/PtValsAlg.cxx,v 1.5 2006/11/11 22:01:59 burnett Exp $
*/


class PtValsAlg : public Algorithm {
public:
    PtValsAlg(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();


private: 
    PointingInfo m_pointing_info;

    StringProperty m_root_tree;
    BooleanProperty m_save_tuple; // set true to save
    StringArrayProperty m_pointingHistory;///< history file name and launch date

    INTupleWriterSvc* m_rootTupleSvc;
    IDataProviderSvc* m_pEventSvc;

    astro::PointingHistory* m_history;

};
//------------------------------------------------------------------------


static const AlgFactory<PtValsAlg>  Factory;
const IAlgFactory& PtValsAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
PtValsAlg::PtValsAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator) 
{
    // declare properties with setProperties calls

    declareProperty("pointing_info_tree_name",  m_root_tree="MeritTuple");
    declareProperty("save_pointing_info",  m_save_tuple=false);
    declareProperty("PointingHistory",   m_pointingHistory); // doublet, filename and launch date

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
    }else if( !m_root_tree.value().empty() ) {
        
        m_pointing_info.setPtTuple(m_rootTupleSvc, m_root_tree.value());
    }

        //set the input file to be used as the pointing database, if used
    if( m_pointingHistory.value().empty()){

        log << MSG::ERROR << "No history file specified" << endreq;
        return StatusCode::FAILURE;

    }
    std::string filename(m_pointingHistory.value()[0]);
    facilities::Util::expandEnvVar(&filename);
    double offset = 0;
    bool horizontalflag(false);
    if( m_pointingHistory.value().size()>1){
        std::string field(m_pointingHistory.value()[1]);
        if(! field.empty() ) { // allow null string
            facilities::Timestamp jt(m_pointingHistory.value()[1]);
            offset = (astro::JulianDate(jt.getJulian())-astro::JulianDate::missionStart())*astro::JulianDate::secondsPerDay;
        }
    }

    if( m_pointingHistory.value().size()>2){
        std::string field(m_pointingHistory.value()[2]);
        horizontalflag =! field.empty();
    }
    log << MSG::INFO << "Loading Pointing History File : " << filename 
        << " with MET offset "<< offset <<  endreq;
    if( horizontalflag){
        log << MSG::INFO << "Will override x-direction to be horizontal"<<endreq;
    }
    m_history = new astro::PointingHistory(filename, offset);

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
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    //
    // Purpose: set tuple items
 
   SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);

   double etime(header->time());
    m_pointing_info.set(etime, (*m_history)(etime));

   
    // put pointing stuff into the root tree
    if( m_rootTupleSvc!=0 && !m_root_tree.value().empty()){
        m_rootTupleSvc->storeRowFlag(this->m_root_tree.value(), m_save_tuple);
    }


    // Here the TDS receives the exposure data
    Event::ExposureCol* exposureDBase = new Event::ExposureCol;
    sc=eventSvc()->registerObject(EventModel::MC::ExposureCol , exposureDBase);
    if(sc.isFailure()) {
        log << MSG::ERROR << EventModel::MC::ExposureCol  
            <<" could not be entered into existing data store" << endreq;
        return sc;
    }
    exposureDBase->push_back(m_pointing_info.forTDS());

    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode PtValsAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    return sc;
}


