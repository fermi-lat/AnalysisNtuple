// Implements ntuple writing algorithm

#include "src/ValBase.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

void ValBase::zeroVals()
{
    mapIter it = m_ntupleMap.begin();
    for (it; it!=m_ntupleMap.end(); it++) {
        *(it->second) = 0.0;
    }
}

StatusCode ValBase::addNtupleValue(std::string varName, INTupleWriterSvc* pSvc, std::string tupleName) 
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
    std::string delim     = "\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    if (varName!="") {
        std::cout  << " Run " << m_run << ", event " << m_event << " -- " ;
    } else {
        std::cout  << " Values of the variables for run " << m_run << ", event " 
        << m_event << ":" << std::endl;
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
        std::cout << delim << it->first << delim << ": " << *(it->second) << " ";
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
    // Recover EventHeader Pointer
    //std::cout << " EventSvc Pointer " << pEventSvc << std::endl;
    SmartDataPtr<Event::EventHeader> pEvent(pEventSvc, EventModel::EventHeader);
    if( pEvent->run()!=m_run || pEvent->event()!=m_event) {
        zeroVals();
        sc = calculate();
        //std::cout << "calc done for this event" << std::endl;
        m_run   = pEvent->run();
        m_event = pEvent->event();
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
    mapIter it = m_ntupleMap.begin();
    
    std::string delim     = "\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    std::cout  << " ValsTool called with unknown name: " << delim << varName << delim << std::endl;
    std::cout << " Known names are: " ;
    int count;
    int length = indent.size();
    it = m_ntupleMap.begin();
    for (it, count=0; it!=m_ntupleMap.end(); it++, count++) {
        length += ((it->first).size() + 2*delim.size() + separator.size());
        if(length>78) {
            std::cout <<std::endl << indent ;
            length = indent.size();
        }
        std::cout << delim << it->first << delim << " ";
    }
    std::cout << std::endl;
}


void ValBase::setEventSvc(IDataProviderSvc* svc)
{
    pEventSvc = svc;
}

StatusCode ValBase::calculate() {
    std::cout << "No specific calc routine defined!" << std::endl;
    return StatusCode::SUCCESS;
}

