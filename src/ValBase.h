/** @file ValBase.h
@brief header file for ValBase.cxx
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ValBase.h,v 1.18 2003/09/30 19:01:20 lsrea Exp $
*/

#ifndef ValBase_h
#define ValBase_h

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IIncidentListener.h"

#include "AnalysisNtuple/IValsTool.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>

class IIncidentSvc;
class IDataProviderSvc;

/** @class ValBase
@brief base class for the XxxValsTools

@author Leon Rochester

*/

namespace {
    // for general use among the ValsTools
    //double M_PI = 3.14159265358979323846;
}

class ValBase : public IValsTool,  public AlgTool,  virtual public IIncidentListener
{
public:
    
    typedef std::pair<std::string, double*> valPair;
    typedef std::vector<valPair*> valMap;
    typedef valMap::iterator mapIter;
    typedef valMap::const_iterator constMapIter;
    
    ValBase(const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    ~ValBase(); 
    /// clear map values
    virtual void zeroVals();
    /// add an item to the map
    virtual void addItem(std::string varName, double* pValue);
    /// do calculation if not already done for this event
    virtual StatusCode doCalcIfNotDone();
    /// get a particular value, using ntuple name default forces calculation
    virtual StatusCode getVal(std::string varName, double& value, int check = 0);
    /// get a particular value, using ntuple name, with calc checking (called by AnaTup)
    virtual StatusCode getValCheck(std::string varName, double& value);
    /// output the list of names
    virtual void announceBadName(std::string varName);
    /// output the names and values, either all (default) or just one;
    virtual StatusCode browse(std::string varName = "");
    /// this is called by the incident service at the beginning of an event
    virtual void handle(const Incident& inc);
    /// callback for visitor
    virtual IValsTool::Visitor::eVisitorRet traverse(IValsTool::Visitor * v,
        const bool checkCalc);
    virtual int getCalcCount() { return m_calcCount;}
    
    /// calculate all values; implemented by each XxxValsTool
    virtual StatusCode calculate();
    
    // common initialization
    virtual StatusCode initialize();
    
protected:
    
    /// map containing ntuple names, and pointers to the ntuple variables
    valMap m_ntupleMap;
    /// pointer to incident service
    IIncidentSvc* m_incSvc;
    /// let ValBase handle the pointer to the data service, everyone uses it
    IDataProviderSvc* m_pEventSvc;
    /// flag to signal new event
    bool m_newEvent;
    /// flag to allow an always-calculate call if 0; if 1 checks and sets m_newEvent
    /// if -1 never skips calculation
    int m_check;

    /// count calls to tools
    static int m_calcCount;
};
#endif