/** @file ValBase.cxx
@brief implements all the methods of the XxxValsTools
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ValBase.cxx,v 1.17 2003/05/07 15:05:20 lsrea Exp $
*/

#include "ValBase.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include <algorithm>
#include <cmath>

int ValBase::m_calcCount = 0;

ValBase::ValBase(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
  : AlgTool( type, name, parent ) { }


StatusCode ValBase::initialize()
{
    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    IDataProviderSvc* eventsvc = 0;

    m_newEvent = true;
    m_check = true;

    m_ntupleMap.clear();

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
    
    //set up listener for IncidentSvc
    incsvc->addListener(this, "BeginEvent", 100);
    return sc;
}

ValBase::~ValBase()
{
    mapIter it = m_ntupleMap.begin();
    for ( ; it!=m_ntupleMap.end(); ++it) {
        delete (*it);
    }
}

void ValBase::zeroVals()
{
    mapIter it = m_ntupleMap.begin();
    for ( ; it!=m_ntupleMap.end(); ++it) {
        *((*it)->second) = 0.0;
    }

}

void ValBase::addItem(std::string varName, double* pValue)
{
    valPair* pair = new valPair(varName, pValue);

    m_ntupleMap.push_back(pair);
}


StatusCode ValBase::browse(std::string varName) 
{
    // browse always triggers a calculation, which doesn't reset the m_newEvent flag
    MsgStream log(msgSvc(), name());

    m_check = false;

    log << MSG::INFO << "ValBase::browse called" << endreq;

    std::string delim     = "\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    if(doCalcIfNotDone().isFailure()) {
        m_check = true;
        return StatusCode::FAILURE;
    }
    m_check = true;
    
    if (varName!="") {
        std::cout  << " Variable " ;
    } else {
        std::cout   << " Values of the variables:" << std::endl << indent;
    }
    int length = indent.size();
    constMapIter it = m_ntupleMap.begin();
    for ( ; it!=m_ntupleMap.end(); ++it) {
        valPair* pair = *it;
        if (varName!="" && varName!=pair->first) continue;
        int valLen = (*(pair->second)==0) ? 1 : 13;
        if (fmod(*(pair->second),1.)==0) valLen = 5;
        int deltaL= (pair->first).size() + 2*delim.size() + separator.size() + valLen + 2;
        length += deltaL;
        if(length>78) {
            std::cout <<std::endl << indent ;
            length = indent.size() + deltaL;
        }
        std::cout  << delim << pair->first << delim << ": " << *(pair->second) << separator;
    }
    std::cout << std::endl;
    return StatusCode::SUCCESS;
}

StatusCode ValBase::doCalcIfNotDone()
{
    StatusCode sc = StatusCode::SUCCESS; 

    MsgStream log(msgSvc(), name());

    // do the calculation if the first time through, or m_check is false

    if(m_newEvent || !m_check) {
        if (!m_pEventSvc)  return StatusCode::FAILURE;
        zeroVals();
        ++m_calcCount;
        sc = calculate();
        //std::cout << "calculation done for this event, total = " << m_calcCount << std::endl;
        // only reset the newEvent flag if we're called with the check flag on.
        if(m_check) m_newEvent = false;
    }
    return sc;
}       

StatusCode ValBase::getValCheck(std::string varName, double& value)
{
    // a simple way to force the check
    return getVal(varName, value, true);
}

StatusCode ValBase::getVal(std::string varName, double& value, bool check)
{
    // optional check flag, if true, do check (default is false)

    StatusCode sc = StatusCode::SUCCESS;

    m_check = check;
    
    constMapIter it = m_ntupleMap.begin();
    for ( ; it!=m_ntupleMap.end(); ++it) {
        if ((*it)->first == varName) break;
    }
    
    if (it==m_ntupleMap.end()) { 
        announceBadName(varName); 
        m_check = true;
        return StatusCode::FAILURE;
    } else {
        if(doCalcIfNotDone().isFailure()) {
            m_check = true;
            return StatusCode::FAILURE;
        }
        value = *(*it)->second;
    }
    m_check = true;
    return sc;
}

void ValBase::announceBadName(std::string varName)  
{
    MsgStream log(msgSvc(), name());

    std::string delim     = "\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    std::cout << " ValsTool called with unknown name: " << delim << varName << delim << std::endl;
    std::cout << " Known names are: " ;

    int length = indent.size();
    int count;

    constMapIter it = m_ntupleMap.begin();
    for (count=0; it!=m_ntupleMap.end(); ++it, ++count) {
        valPair* pair = *it;
        length += ((pair->first).size() + 2*delim.size() + separator.size());
        if(length>78) {
            std::cout << std::endl << indent ;
            length = indent.size();
        }
        std::cout << delim << pair->first << delim << " ";
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
        m_calcCount = 0;
    }
}

IValsTool::Visitor::eVisitorRet ValBase::traverse(IValsTool::Visitor* v,
                                                  const bool checkCalc)
{
    IValsTool::Visitor::eVisitorRet ret = IValsTool::Visitor::DONE;

    if (checkCalc) {
        if(doCalcIfNotDone().isFailure()) return IValsTool::Visitor::ERROR;
    }

    constMapIter it = m_ntupleMap.begin();
    for ( ; it!=m_ntupleMap.end(); ++it) {
        valPair* pair = *it;
        ret = v->analysisValue(pair->first, *(pair->second));
        if (ret!= IValsTool::Visitor::CONT) return ret;
    }
    return IValsTool::Visitor::DONE;
}
