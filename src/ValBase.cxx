// Implements ntuple writing algorithm

#include "ValBase.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"

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

StatusCode ValBase::fillNtuple(INTupleWriterSvc* pSvc, 
                               std::string tupleName,
                               std::string varName) 
{
    StatusCode sc = StatusCode::SUCCESS;
    
    mapIter it = m_ntupleMap.find(varName);
    if (it==m_ntupleMap.end()) {
        announceBadName(varName);
        return StatusCode::FAILURE;
    } else {
        sc = doCalcIfNotDone();
        double value = *(it->second);
        if ((sc = pSvc->addItem(tupleName.c_str(), varName.c_str(), value)).isFailure()) return sc;
    }
    return sc;
}

void ValBase::browseValues(std::string varName) 
{
    MsgStream log(msgSvc(), name());

    std::string delim     = "\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    if (varName!="") {
        std::cout  << " Variable " ;
    } else {
        std::cout   << " Values of the variables:" << std::endl << indent;
    }
    int length = indent.size();
    mapIter it = m_ntupleMap.begin();
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
}

StatusCode ValBase::fillNtuple (INTupleWriterSvc* pSvc, std::string tupleName) 
{
    StatusCode sc = StatusCode::SUCCESS;
    
    sc = doCalcIfNotDone();
    mapIter it = m_ntupleMap.begin();
    for (it; it!=m_ntupleMap.end(); it++) {
        std::string varName = it->first;
        double value       = *(it->second);
        if ((sc = pSvc->addItem(tupleName.c_str(), varName.c_str(), value)).isFailure()) return sc;
    }
    return sc;
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
        std::cout << "calculation done for this event" << std::endl;
        m_newEvent = false;
    }
    return sc;
}       


StatusCode ValBase::getVal(std::string varName, double& value) {
    
    StatusCode sc = StatusCode::SUCCESS;
    
    mapIter it = m_ntupleMap.find(varName);
    
    if (it==m_ntupleMap.end()) { 
        announceBadName(varName); 
        return StatusCode::FAILURE;
    } else {
        sc = doCalcIfNotDone();
        value = *(it->second);
    }
    return sc;
}

void ValBase::announceBadName(std::string varName) {

    MsgStream log(msgSvc(), name());

    mapIter it = m_ntupleMap.begin();
    
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
        std::cout << "handle called at BeginEvent" << std::endl;
        m_newEvent = true;
    }
}

