// Implements ntuple writing algorithm

#include "ValBase.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

ValBase::ValBase(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
  : AlgTool( type, name, parent )
  , m_newEvent(true)
  , m_handleSet(false)
{ 
    m_ntupleMap.clear();
 
}


StatusCode ValBase::initialize()
{
    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    IDataProviderSvc* eventsvc = 0;

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::FAILURE;

    log << MSG::INFO << "ValBase is initializing" << endreq;

    if (serviceLocator()) {
        sc = serviceLocator()->service( "IncidentSvc", incsvc, true );
        if(sc.isFailure()){
            log << MSG::ERROR << "Could not find IncidentSvc" << std::endl;
            return sc;
        }
        m_incSvc = incsvc;
        
        sc = serviceLocator()->service( "EventDataSvc", eventsvc, true );
        if(sc.isFailure()){
            log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
            return sc;
        }
        m_pEventSvc = eventsvc;
    }
    
    //incsvc->addListener(this, "BeginEvent", 100);
    return sc;
}



void ValBase::zeroVals()
{
    mapIter it = m_ntupleMap.begin();
    for (it; it!=m_ntupleMap.end(); it++) {
        *(it->second) = 0.0;
    }
}

void ValBase::addItem(std::string varName, double* pValue)
{
    m_ntupleMap[varName] = pValue;
}

StatusCode ValBase::browse(std::string varName) 
{
    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "ValBase::browse called" << endreq;

    std::string delim     = "\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    if(doCalcIfNotDone().isFailure()) return StatusCode::FAILURE;
    
    if (varName!="") {
        std::cout  << " Variable " ;
    } else {
        std::cout   << " Values of the variables:" << std::endl << indent;
    }
    int length = indent.size();
    constMapIter it = m_ntupleMap.begin();
    for (it; it!=m_ntupleMap.end(); it++) {
        if (varName!="" && varName!=it->first) continue;
        length += (it->first).size() + 2*delim.size() + separator.size() + 15;
        if(length>78) {
            std::cout <<std::endl << indent ;
            length = indent.size();
        }
        std::cout  << delim << it->first << delim << ": " << *(it->second) << " ";
    }
    std::cout << std::endl;
    return StatusCode::SUCCESS;
}

StatusCode ValBase::doCalcIfNotDone()
{
    StatusCode sc = StatusCode::SUCCESS; 

    MsgStream log(msgSvc(), name());
   
#if 1  
    // kludge to get around multiple initializations of tool
    // too many calls to addListener otherwise!
    if(!m_handleSet) {
        m_incSvc->addListener(this, "BeginEvent", 100);
        m_handleSet = true;
    }
#endif

    if(m_newEvent) {
        zeroVals();
        sc = calculate();
        //std::cout << "calculation done for this event" << std::endl;
        m_newEvent = false;
    }
    return sc;
}       


StatusCode ValBase::getVal(std::string varName, double& value)
{
    StatusCode sc = StatusCode::SUCCESS;
    
    constMapIter it = m_ntupleMap.find(varName);
    
    if (it==m_ntupleMap.end()) { 
        announceBadName(varName); 
        return StatusCode::FAILURE;
    } else {
        if(doCalcIfNotDone().isFailure()) return StatusCode::FAILURE;
        value = *(it->second);
    }
    return sc;
}

void ValBase::announceBadName(std::string varName)  
{
    MsgStream log(msgSvc(), name());

    mapIter it;
    
    std::string delim     = "\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    std::cout << " ValsTool called with unknown name: " << delim << varName << delim << std::endl;
    std::cout << " Known names are: " ;
    int count;
    int length = indent.size();
    it = m_ntupleMap.begin();
    for (it, count=0; it!=m_ntupleMap.end(); it++, count++) {
        length += ((it->first).size() + 2*delim.size() + separator.size());
        if(length>78) {
            std::cout << std::endl << indent ;
            length = indent.size();
        }
        std::cout << delim << it->first << delim << " ";
    }
    std::cout << std::endl;
}

StatusCode ValBase::calculate() {

    MsgStream log(msgSvc(), name());

    std::cout << "No specific calc routine defined!" << std::endl;
    return StatusCode::SUCCESS;
}

void ValBase::handle(const Incident & inc) 
{    
    MsgStream log(msgSvc(), name());

    if(inc.type()=="BeginEvent") {
        //std::cout << "handle called at BeginEvent" << std::endl;
        m_newEvent = true;
    }
}

ValsVisitor::eVisitorRet ValBase::traverse(ValsVisitor* v)
{
    StatusCode sc;

    ValsVisitor::eVisitorRet ret = ValsVisitor::DONE;
    constMapIter it = m_ntupleMap.begin();

    if(doCalcIfNotDone().isFailure()) return ValsVisitor::ERROR;
   
    for (it; it!=m_ntupleMap.end(); it++) {
        double value = *(it->second);
        ret = v->analysisValue(it->first, value);
        if (ret!= ValsVisitor::CONT) return ret;
    }
    
    return ValsVisitor::DONE;
}




