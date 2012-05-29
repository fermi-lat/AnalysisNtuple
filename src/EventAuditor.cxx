
#include "GaudiKernel/Auditor.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IChronoStatSvc.h"
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
* $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/EventAuditor.cxx,v 1.8 2011/12/12 20:36:35 heather Exp $
*/

class EventAuditor : virtual public Auditor
{
public:

    EventAuditor( const std::string &, ISvcLocator * ) ;
    virtual StatusCode initialize() ;
    virtual StatusCode finalize() ;
    virtual ~EventAuditor() ;

    virtual void before(StandardEventType, INamedInterface*);
    virtual void after(StandardEventType, INamedInterface*, const StatusCode&);

private :

    INTupleWriterSvc* m_rootTupleSvc;
    StringProperty m_root_tree;
    bool m_save_tuple;


    //! names of algorithms to audit
    std::vector<std::string>  m_algoNames ;
    int m_nAlgs;

    // some pointers to services  
    MsgStream * m_log ;
    IChronoStatSvc * m_chronoSvc ;

    // internals
    std::map<IAlgorithm*,      bool> m_algoStatus ;
    std::map<INamedInterface*, bool> m_nameStatus ;

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
<tr><td> AudCalCluster
<td>F<td>   Time for the calorimeter clustering, part of AudCal1
<tr><td> AudCal1
<td>F<td>   Time for the calorimeter part of the reconstruction (pass 1)
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
    algoNamesVec.push_back("Cal1"); 
    algoNamesVec.push_back("Tkr");
    algoNamesVec.push_back("Cal2");
    algoNamesVec.push_back("CalCluster"); 
    
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

void EventAuditor::before(StandardEventType type, INamedInterface* namedInterface)
{
    // Only running on "execute"
    if (type != IAuditor::Execute) 
    {
        return;
    }

    // upgrade m_algoStatus
    std::string thisName = namedInterface->name();
    if(thisName=="Event") {
        // initialize the timers
        m_timeVals.assign(m_nAlgs, -1.0);
    } 

    std::map<INamedInterface*, bool>::iterator itr ;
    itr = m_nameStatus.find(namedInterface);
    if ( itr == m_nameStatus.end() )
    {
        const std::vector< std::string > & algoNames = m_algoNames ;
        std::vector< std::string >::const_iterator algoName ;

        algoName = std::find(algoNames.begin(), algoNames.end(), thisName);
        if(algoName==algoNames.end()) {
            m_nameStatus[namedInterface] = false ; 
            return;
        } else {
            m_nameStatus[namedInterface] = true ;
        }
    }

    if(m_nameStatus[namedInterface]) {
        m_chronoSvc->chronoStart("EA_"+thisName) ;
        (*m_log) << MSG::DEBUG << "start " << thisName << std::endl;
    }
    return;
}

void EventAuditor::after(StandardEventType type, INamedInterface* namedInterface, const StatusCode& sc)
{
    std::map<INamedInterface*,bool>::iterator itr = m_nameStatus.find(namedInterface);

    // look if the algo is under monitoring
    if (itr == m_nameStatus.end() || itr->second == false || type != IAuditor::Execute)
    { 
        return;
    }

    std::string thisName = namedInterface->name();

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
    return;
}


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
