/** @file ValBase.cxx
@brief implements all the methods of the XxxValsTools
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ValBase.cxx,v 1.32 2007/03/14 22:12:21 lsrea Exp $
*/

#include "ValBase.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include <algorithm>
#include <cassert>

const int ValBase::s_badVal = -9999;

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
      m_check = CHECK;
      m_calcCount = 0;

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
        delete (*it)->second;
        delete (*it);
    }
}

void ValBase::zeroVals()
{
    mapIter it = m_ntupleMap.begin();
    for ( ; it!=m_ntupleMap.end(); ++it) {
        TypedPointer* ptr = ((*it)->second);
        void* vPtr = ptr->getPointer();
        valType type = ptr->getType();
        int dim = ptr->getDim();

        int i;
        for (i=0; i<dim; ++i) {
            switch (type) {
            case FLOAT: 
                *(reinterpret_cast<float*>(vPtr)+i) = 0.0;
                //std::cout << "After zero (F): " << *(reinterpret_cast<float*>(vPtr)) << std::endl;
                break;
            case DOUBLE:
                *(reinterpret_cast<double*>(vPtr)+i) = 0.0;
                break;
            case INT:
                *(reinterpret_cast<int*>(vPtr)+i) = 0;
                break;
            case UINT:
                *(reinterpret_cast<unsigned int*>(vPtr)+i) = 0;
                //std::cout << "After zero (I): " << *(reinterpret_cast<int*>(vPtr)) << std::endl;
                break;
            case STRING:
                if (i==0) strcpy(reinterpret_cast<char*>(vPtr), "_");
                //std::cout << "After zero (S): *" << reinterpret_cast<char*>(vPtr) 
                //    << "*" << std::endl;
                break;
            }
        }
    }
}

std::string ValBase::getFullName(std::string varName, int dim)
{
    char buffer[6];
    sprintf(buffer, "[%i]", dim);
    std::string fullName = varName+buffer;
    return fullName;
}

bool ValBase::getArrayArg(std::string varName, std::string& baseName, int& arg)
{
    MsgStream log(msgSvc(), name());
    
    arg = -1;
    baseName = varName;
    bool hasArg = false;

    int pos1, pos2;
    pos1 = varName.find("[");
    if (pos1!=-1) {
        pos2 = varName.find("]");
        if (pos2<pos1+2) {
            log << MSG::ERROR << "variable " << varName << " out of range or malformed" << endreq;
            log << MSG::ERROR << "Character positions " << pos1 << " " << pos2 << endreq;
            assert(pos2<pos1+2);
            return hasArg;
        } else {
            std::string strDim = varName.substr(pos1+1, pos2-pos1-1);
            arg = atoi(strDim.c_str());
            if(arg<0) {
                log << MSG::ERROR << "variable " << varName << " out of range or malformed" << endreq;
                log << MSG::ERROR << "argString = " << strDim << ", arg = " << arg << endreq;
                assert(arg<0);
                return hasArg;
            }
            baseName = varName.substr(0,pos1);
            hasArg = true;
            return hasArg;
        }
    } else {
        arg = 1;
        return hasArg;
    }
}

void ValBase::addItem(std::string varName, double* pValue)
{
    std::string baseName;
    int dim;
    /*bool hasArg =*/ getArrayArg(varName, baseName, dim);
    TypedPointer* ptr = new TypedPointer(DOUBLE, (void*) pValue, dim);
    valPair* pair = new valPair(baseName, ptr);

    m_ntupleMap.push_back(pair);
}

void ValBase::addItem(std::string varName, float* pValue)
{
    std::string baseName;
    int dim;
    /*bool hasArg =*/ getArrayArg(varName, baseName, dim);
    TypedPointer* ptr = new TypedPointer(FLOAT, (void*) pValue, dim);
    valPair* pair = new valPair(baseName, ptr);

    m_ntupleMap.push_back(pair);
}
void ValBase::addItem(std::string varName, int* pValue)
{
    std::string baseName;
    int dim;
    /*bool hasArg =*/ getArrayArg(varName, baseName, dim);
    TypedPointer* ptr = new TypedPointer(INT, (void*) pValue, dim);
    valPair* pair = new valPair(baseName, ptr);

    m_ntupleMap.push_back(pair);
}

void ValBase::addItem(std::string varName, unsigned int* pValue)
{
    std::string baseName;
    int dim;
    /*bool hasArg =*/ getArrayArg(varName, baseName, dim);
    TypedPointer* ptr = new TypedPointer(UINT, (void*) pValue, dim);
    valPair* pair = new valPair(baseName, ptr);

    m_ntupleMap.push_back(pair);
}

void ValBase::addItem(std::string varName, char* pValue)
{
    std::string baseName;
    int dim;
    bool hasArg = getArrayArg(varName, baseName, dim);
    if(hasArg) {
        dim = 1;
        std::cout << std::endl << "AnalysisNtuple/varBase " << std::endl;
        std::cout << "     " << varName << ": arrays of strings not allowed" << std::endl;
        std::cout << "     first one will be used... " << std::endl << std::endl;
    }

    TypedPointer* ptr = new TypedPointer(STRING, (void*) pValue, dim);
    valPair* pair = new valPair(baseName, ptr);

    m_ntupleMap.push_back(pair);
}

StatusCode ValBase::browse(MsgStream log, std::string varName0) 
{
    // browse always triggers a calculation, which doesn't reset the m_newEvent flag

    m_check = CALC;

    //log << MSG::INFO << "ValBase::browse called" << endreq;

    std::string delim     = ""; //"\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    if(doCalcIfNotDone().isFailure()) {
        m_check = CHECK;
        return StatusCode::FAILURE;
    }
    m_check = CHECK;

    std::string varName;
    int element;
    bool hasArg = getArrayArg(varName0, varName, element);
    if (!hasArg) element = 0;
    bool doAll = !hasArg;

    
    if (varName!="") {
        log << MSG::INFO << " Variable " ;
    } else {
        log   << MSG::INFO << " Values of the variables:" << endreq << indent;
    }
    int length = indent.size();
    constMapIter it = m_ntupleMap.begin();
    
    for ( ; it!=m_ntupleMap.end(); ++it) {
        valPair* pair = *it;
        if (varName!="" && varName!=pair->first) continue;
        int valLen = 13;
        int deltaL= (pair->first).size() + 2*delim.size() + separator.size() + valLen + 2;
        length += deltaL;
        if(length>78) {
            log << MSG::INFO << endreq << indent ;
            length = indent.size() + deltaL;
        }
        log << delim << pair->first << delim << ": " ;

        TypedPointer* ptr = (*it)->second;
        valType type = ptr->getType();
        int dim = ptr->getDim();
        if(element>=dim || element<0) {
            log << MSG::ERROR << "Browse: error in arg: " << element << endreq;
            assert(element>=dim || element<0);
        }

        int start = (doAll ? 0 : element);
        int end   = (doAll ? dim-1 : element);
        void* vPtr = ptr->getPointer();

        int i;
        if (doAll && dim>1) log << "(";
        
        for (i=start; i<=end; ++i) {

            switch (type) {
            case FLOAT: 
                log << *(reinterpret_cast<float*>(vPtr)+i);
                break;
            case DOUBLE:
                log << *(reinterpret_cast<double*>(vPtr)+i);
                break;
            case INT:
                log << *(reinterpret_cast<int*>(vPtr)+i);
                break;
            case UINT: 
                log << *(reinterpret_cast<unsigned int*>(vPtr)+i);
                break;
            default:
                break;
            }

            if(doAll && dim>1) {
                log << (i==dim-1 ? ")" : ",");
            }
            log << separator;
        }
        
    }
    log << endreq;
    return StatusCode::SUCCESS;
}

StatusCode ValBase::doCalcIfNotDone()
{
    StatusCode sc = StatusCode::SUCCESS; 

    MsgStream log(msgSvc(), name());

    // if NOCALC means don't do the calculation
    // if CALC   means always do the calculation
    // if CHECK  means do the calculation if not already done, and reset m_newEvent

    if(m_check!=NOCALC) {
        if(m_newEvent || m_check==CALC) {
            if (!m_pEventSvc)  return StatusCode::FAILURE;
            zeroVals();
            ++m_calcCount;
            sc = calculate();
            //std::cout  << m_calcCount << " calculations so far" << std::endl;
            // only reset the newEvent flag if we're called with the check flag on.
            if(m_check==CHECK) m_newEvent = false;
        }
    }
    return sc;
}       

StatusCode ValBase::getValCheck(std::string varName, double& value)
{
    // a simple way to force the check
    return getVal(varName, value, CHECK);
}
StatusCode ValBase::getValCheck(std::string varName, float& value)
{
    // a simple way to force the check
    return getVal(varName, value, CHECK);
}
StatusCode ValBase::getValCheck(std::string varName, int& value)
{
    // a simple way to force the check
    return getVal(varName, value, CHECK);
}

StatusCode ValBase::getValCheck(std::string varName, unsigned int& value)
{
    // a simple way to force the check
    return getVal(varName, value, CHECK);
}

StatusCode ValBase::getValCheck(std::string varName, std::string& value)
{
    // a simple way to force the check
    return getVal(varName, value, CHECK);
}

StatusCode ValBase::getTypedPointer(std::string varName, TypedPointer*& ptr, int check)
{
    // optional check flag

    StatusCode sc = StatusCode::SUCCESS;

    m_check = check;
    
    constMapIter it = m_ntupleMap.begin();
    for ( ; it!=m_ntupleMap.end(); ++it) {
        if ((*it)->first == varName) break;
    }
    
    if (it==m_ntupleMap.end()) { 
        announceBadName(varName); 
        m_check = CHECK;
        return StatusCode::FAILURE;
    } else {
        if(doCalcIfNotDone().isFailure()) {
            m_check = CHECK;
            return StatusCode::FAILURE;
        }
        ptr = (*it)->second;
    }
    m_check = CHECK;
    return sc;
}

StatusCode ValBase::getVal(std::string varName, std::string& value, int check)
{
    MsgStream log(msgSvc(), name());

    char buffer[80];
    TypedPointer* ptr = 0;
    value = "";

    std::string baseName;
    int element;
    bool hasArg = getArrayArg(varName, baseName, element);
    if(!hasArg) element = 0;
    StatusCode sc = getTypedPointer(baseName, ptr, check);
    if(sc.isFailure()) return sc;

    void* vPtr = ptr->getPointer();

    int dim = ptr->getDim();
    if (element>=dim || element<0) 
    {
        log << MSG::ERROR << "GetVal: element " << varName << " out of range" << endreq;
            assert(element<dim && element>-1);
    }

    valType type = ptr->getType();
    
    if(sc.isSuccess()) {
        if (type==FLOAT) { sprintf(buffer, "%g", 
            *(reinterpret_cast<float*>(vPtr)+element));}
        else if (type==DOUBLE)  { sprintf(buffer, "%g", 
            *(reinterpret_cast<double*>(vPtr)+element));}
        else if (type==INT)     { sprintf(buffer, "%i", 
            *(reinterpret_cast<int*>(vPtr)+element));}
        else if (type==UINT)    { sprintf(buffer, "%u", 
            *(reinterpret_cast<unsigned int*>(vPtr)+element));}
    }
    value = std::string(buffer);
    return sc;
}

StatusCode ValBase::getVal(std::string varName, int& value, int check)
{
    MsgStream log(msgSvc(), name());
    TypedPointer* ptr = 0;
    value = 0;

    std::string baseName;
    int element;
    bool hasArg = getArrayArg(varName, baseName, element);
    if(!hasArg) element = 0;

    StatusCode sc = getTypedPointer(baseName, ptr, check);
    if(sc.isFailure()) return sc;

    int dim = ptr->getDim();
    if (element>=dim || element<0) 
    {
        log << MSG::ERROR << "GetVal: element " << varName << " out of range" << endreq;
            assert(element<dim && element>-1);
    }

    value = *(reinterpret_cast<int*>(ptr->getPointer())+element);

    return sc;
}

StatusCode ValBase::getVal(std::string varName, unsigned int& value, int check)
{
    MsgStream log(msgSvc(), name());
    TypedPointer* ptr = 0;
    value = 0;

    std::string baseName;
    int element;
    bool hasArg = getArrayArg(varName, baseName, element);

    StatusCode sc = getTypedPointer(baseName, ptr, check);
    if(sc.isFailure()) return sc;

    if(!hasArg) element = 0;
    int dim = ptr->getDim();
    if (element>=dim || element<0) 
    {
        log << MSG::ERROR << "GetVal: element " << varName << " out of range" << endreq;
            assert(element<dim && element>-1);
    }

    value = *(reinterpret_cast<unsigned int*>(ptr->getPointer())+element);
   
    return sc;
}
StatusCode ValBase::getVal(std::string varName, double& value, int check)
{
    MsgStream log(msgSvc(), name());
    TypedPointer* ptr = 0;
    value = 0.0;

    std::string baseName;
    int element;
    bool hasArg = getArrayArg(varName, baseName, element);
    if(!hasArg) element = 0;

    StatusCode sc = getTypedPointer(baseName, ptr, check);
    if(sc.isFailure()) return sc;
    int dim = ptr->getDim();
    if (element>=dim || element<0) 
    {
        log << MSG::ERROR << "GetVal: element " << varName << " out of range" << endreq;
            assert(element<dim && element>-1);
    }

    value = *(reinterpret_cast<double*>(ptr->getPointer())+element);

    return sc;
}

StatusCode ValBase::getVal(std::string varName, float& value, int check)
{
    MsgStream log(msgSvc(), name());
    TypedPointer* ptr = 0;
    value = 0.0;

    std::string baseName;
    int element;
    bool hasArg = getArrayArg(varName, baseName, element);

    StatusCode sc = getTypedPointer(baseName, ptr, check);
    if(sc.isFailure()) return sc;

    int dim = ptr->getDim();
    if(!hasArg) element = 0;
    if (element>=dim || element<0) 
    {
        log << MSG::ERROR << "GetVal: element " << varName << " out of range" << endreq;
            assert(element<dim && element>-1);
    }

    value = *(reinterpret_cast<float*>(ptr->getPointer())+element);

    return sc;
}

void ValBase::announceBadName(std::string varName)  
{
    MsgStream log(msgSvc(), name());

    std::string delim     = ""; //"\"";
    std::string separator = " ";
    std::string indent    = "    ";
    
    std::cout << " ValsTool called with unknown name: " << delim << varName << delim << std::endl;
    /*
    std::cout << " Known names are: " ;

    int length = indent.size();
    int count;

    constMapIter it = m_ntupleMap.begin();
    for (count=0; it!=m_ntupleMap.end(); ++it, ++count) {
        valPair* pair = *it;
        std::string thisName = pair->first;
        int dim = (pair->second)->getDim();
        char result[80];
        sprintf(result, "%i",dim);
        if (dim>1) thisName = thisName+"["+result+"]";
        length += (thisName.size() + 2*delim.size() + separator.size());
        if(length>78) {
            std::cout << std::endl << indent ;
            length = indent.size();
        }
        std::cout << delim << thisName << delim << " ";
    }
    std::cout << std::endl;
    */
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
        TypedPointer* ptr = pair->second;
        valType type = ptr->getType();
        void* vPtr = ptr->getPointer();
        std::string varName = pair->first;

        // here's where we need to construct the varName from the baseName and the dim
        std::string fullName = varName;
        int dim = ptr->getDim();
        if (dim>1) {
            fullName = getFullName(varName, dim);
        }

        switch (type) {
            case FLOAT:
                ret = v->analysisValue(fullName, *(reinterpret_cast<float*>(vPtr)));
                //std::cout << "analysisValue returns: " << *(reinterpret_cast<float*>(vPtr)) << std::endl;
                break;
            case DOUBLE:
                ret = v->analysisValue(fullName, *(reinterpret_cast<double*>(vPtr)));
                break;
            case UINT:
                ret = v->analysisValue(fullName, *(reinterpret_cast<unsigned int*>(vPtr)));
                break;
            case INT:
                ret = v->analysisValue(fullName, *(reinterpret_cast<int*>(vPtr)));
                //std::cout << "analysisValue returns: " << *(reinterpret_cast<int*>(vPtr)) << std::endl;
               break;
            case STRING:
                ret = v->analysisValue(fullName, (reinterpret_cast<char*>(vPtr)));
                //std::cout << "analysisValue returns: " << reinterpret_cast<char*>(vPtr) << std::endl;
                break;
        }

        if (ret!= IValsTool::Visitor::CONT) return ret;
    }
    return IValsTool::Visitor::DONE;
}
