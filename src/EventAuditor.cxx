
#include "GaudiKernel/Auditor.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Chrono.h"
#include "GaudiKernel/IAlgorithm.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <map>
#include <algorithm>

/**   
* @class EventAuditor
*
* Experiment the monitoring of algorithms.
*
* $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/EventAuditor.cxx,v 1.5 2010/01/01 02:06:13 lsrea Exp $
*/

class EventAuditor : virtual public Auditor
{
public:

    EventAuditor( const std::string &, ISvcLocator * ) ;
    virtual StatusCode initialize() ;
    virtual StatusCode finalize() ;
    virtual ~EventAuditor() ;

    //virtual StatusCode beforeInitialize( IAlgorithm * ) ;
    //virtual StatusCode afterInitialize( IAlgorithm * ) ;
    virtual StatusCode beforeExecute( IAlgorithm * ) ;
    virtual StatusCode afterExecute( IAlgorithm * ) ;
    //virtual StatusCode beforeFinalize( IAlgorithm * ) ;
    //virtual StatusCode afterFinalize( IAlgorithm * ) ;

private :

    INTupleWriterSvc* m_rootTupleSvc;
    StringProperty m_root_tree;
    bool m_save_tuple;


    //! names of algorithms to audit
    //StringArrayProperty  m_algoNames ;
    std::vector<std::string>  m_algoNames ;
    int m_nAlgs;

    // some pointers to services  
    MsgStream * m_log ;
    IChronoStatSvc * m_chronoSvc ;

    // internals
    std::map< IAlgorithm *, bool > m_algoStatus ;

    std::vector<float> m_timeVals;
} ;

//==========================================================
// construction
//==========================================================

DECLARE_AUDITOR_FACTORY(EventAuditor) ;

/** @page anatup_vars 

@section timing  EventAuditor Times for Selected Algorithms/Sequences

Times for selected entities are given in seconds. The granularity is ~16msec (60Hz). The minimum time is recorded if the clock changes between the start and stop calls of the process.

Note: the time for the first event in a run includes initialization overhead.

The variables are named "Aud"+Entity Name. The current default set is:

<table>
<tr><th> Variable <th> Type <th> Description                                        
<tr><td> AudEvent 
<td>F<td>   Time for the entire event
<tr><td> AudGeneration 
<td>F<td>   Time for the generation
<tr><td> AudReconstruction
<td>F<td>   Time for the reconstruction
<tr><td> AudTkr 
<td>F<td>   Time for the tracker part of the reconstruction (pass 1)
<tr><td> AudCal2
<td>F<td>   Time for the calorimeter part of the reconstruction (pass 2)
</table>

*/


EventAuditor::EventAuditor( const std::string & name, ISvcLocator * svcLocator )
: Auditor(name,svcLocator)
{
    std::vector<std::string> algoNamesVec ;
    algoNamesVec.push_back("Event");
    algoNamesVec.push_back("Generation");
    algoNamesVec.push_back("Reconstruction");
    algoNamesVec.push_back("Tkr");
    algoNamesVec.push_back("Cal2");

    declareProperty("algoNames",m_algoNames=algoNamesVec);
    declareProperty("tree_name",  m_root_tree="MeritTuple");
}

StatusCode EventAuditor::initialize()
{
    Auditor::initialize() ;
    StatusCode sc = StatusCode::SUCCESS;

    setProperties();
    std::vector<std::string>& algoNames = m_algoNames;
    // make sure that Event is one of the monitored algs
    // this is so that we know when there's a new event
    // and can set the times to a default value (-1 sec)
    std::vector< std::string >::const_iterator algoName ;
    algoName = std::find(algoNames.begin(), algoNames.end(), "Event");
    if(algoName==algoNames.end()) {
        algoNames.insert(algoNames.begin(),"Event");
    }

    m_nAlgs = algoNames.size();
    m_timeVals.assign(m_nAlgs, -1.0);

    // message stream
    IMessageSvc * messageSvc ;
    service("MessageSvc",messageSvc,true) ;
    m_log = new MsgStream(messageSvc,name()) ;

    // chrono stat svc
    if (service("ChronoStatSvc",m_chronoSvc,true).isFailure())
    {
        (*m_log)<<MSG::ERROR<<"Could not find TkrReconSvc"<<endreq ;
        return StatusCode::FAILURE ;
    }

    // get a pointer to RootTupleSvc
    if( (service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
        (*m_log) << MSG::ERROR << " RootTupleSvc is not available" << endreq;
        m_rootTupleSvc=0;
        return StatusCode::FAILURE;
    }

    if( m_rootTupleSvc==0 ) return sc;

    std::string tname = m_root_tree.value();

    int i;
    for(i=0;i<m_nAlgs;++i) {
        m_rootTupleSvc->addItem(tname, "Aud"+m_algoNames[i] , &m_timeVals[i]);
    }

    (*m_log)<<MSG::DEBUG<<"initialize() "<<endreq ;
    return StatusCode::SUCCESS ;
}

//==========================================================
// events loop
//==========================================================

//StatusCode EventAuditor::beforeInitialize( IAlgorithm * alg )
//{
//    (*m_log)<<MSG::DEBUG<<"beforeInitialize("<<alg->name()<<") "<<endreq ;
//    return StatusCode::SUCCESS ;
//}
//
//
//StatusCode EventAuditor::afterInitialize( IAlgorithm * algo )
//{
//    (*m_log)<<MSG::DEBUG<<"afterInitialize("<<algo->name()<<") "<<endreq ;
//    return StatusCode::SUCCESS ;
//}

StatusCode EventAuditor::beforeExecute( IAlgorithm * algo )
{
    // upgrade m_algoStatus
    std::string thisName = algo->name();
    if(thisName=="Event") {
        // initialize the timers
        m_timeVals.assign(m_nAlgs, -1.0);
    } 

    std::map< IAlgorithm *, bool >::iterator itr ;
    itr = m_algoStatus.find(algo) ;
    if ( itr == m_algoStatus.end() )
    {
        const std::vector< std::string > & algoNames = m_algoNames ;
        std::vector< std::string >::const_iterator algoName ;

        algoName = std::find(algoNames.begin(), algoNames.end(), thisName);
        if(algoName==algoNames.end()) {
            m_algoStatus[algo] = false ; 
            return StatusCode::SUCCESS ;
        } else {
            m_algoStatus[algo] = true ;
        }
    }

    if(m_algoStatus[algo]) {
        m_chronoSvc->chronoStart("EA_"+thisName) ;
        (*m_log) << MSG::DEBUG << "start " << thisName << std::endl;
    }
    return StatusCode::SUCCESS ;
}

StatusCode EventAuditor::afterExecute( IAlgorithm * algo )
{
    // look if the algo is under monitoring
    if (m_algoStatus[algo]==false)
    { return StatusCode::SUCCESS ; }

    std::string thisName = algo->name();

    // stop chrono
    m_chronoSvc->chronoStop("EA_"+thisName) ;

    // retrieve and log the last time interval
    
    IChronoStatSvc::ChronoTime delta
        = m_chronoSvc->chronoDelta("EA_"+thisName,IChronoStatSvc::USER) ;

    float fDelta = static_cast<float>(delta)*0.000001;
    (*m_log) << MSG::DEBUG << thisName <<" user time: " << fDelta << " sec" << endreq ;

    const std::vector< std::string > & algoNames = m_algoNames ;
    int i;
    for(i=0;i<m_nAlgs;++i) {
        if(thisName==algoNames[i]) {
            m_timeVals[i] = fDelta;
            break;
        }
    }

/* HMK Apr172008
    Do not want EventAuditor forcing events to be stored to the MeritTuple
    without meeting our trigger requirements.  Other algorithms will set
    this flag for us, such as AnalysisNtupleAlg

    m_save_tuple = true;
    if( m_rootTupleSvc!=0 && !m_root_tree.value().empty()){
        m_rootTupleSvc->storeRowFlag(m_root_tree.value(), m_save_tuple);
    }
*/

    return StatusCode::SUCCESS ;
}

//StatusCode EventAuditor::beforeFinalize( IAlgorithm * alg )
// {
//     (*m_log) << MSG::DEBUG << "beforeFinalize("<<alg->name()<<") "<< endreq ;
//  return StatusCode::SUCCESS ;
// }
//
//
//StatusCode EventAuditor::afterFinalize( IAlgorithm * alg )
// {
//     (*m_log) << MSG::DEBUG << "afterFinalize("<<alg->name()<<") "<< endreq ;
//  return StatusCode::SUCCESS ;
// }
//

//==========================================================
// destruction
//==========================================================

StatusCode EventAuditor::finalize()
{
    (*m_log) << MSG::DEBUG <<"finalize() "<< endreq ;
//    return Auditor::finalize() ;
  return StatusCode::SUCCESS;

}

EventAuditor::~EventAuditor()
{}



