#ifndef ValBase_h
#define ValBase_h

#include "GaudiKernel/AlgTool.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IIncidentListener.h"

#include "AnalysisNtuple/IValsTool.h"
#include <string>
#include <map>

/** @class ValBase
  @brief 

  */
class ValBase : public IValsTool,  public AlgTool,  virtual public IIncidentListener
{
public:

typedef std::map<std::string, double*> valMap;
typedef valMap::iterator mapIter;

    ValBase(const std::string& type, 
            const std::string& name, 
            const IInterface* parent);
    
    ~ValBase() {return;}
    
    /// clear map values
    virtual void zeroVals();
    /// fill ntuple values into tupleSvc
    virtual StatusCode fillNtuple(INTupleWriterSvc* pSvc, std::string tupleName);
    /// put one value into the ntuple
    virtual StatusCode fillNtuple(INTupleWriterSvc* pSvc,
        std::string tupleName, std::string varName);
    /// check if calculation is already done for this event
    virtual StatusCode doCalcIfNotDone();
    /// get a particular value, using ntuple name
    virtual StatusCode getVal(std::string varName, double& value);
    /// output the list of names
    virtual void announceBadName(std::string varName);
    /// output the names and values, either all (default) or just one;
    virtual void browseValues(std::string varName = "");
    /// this is called by the incident service at the beginning of an event
    virtual void handle(const Incident& inc);
    
    /// calculate all values, over-ridden by XxxValsTool
    virtual StatusCode calculate();

    // common initialization for subclasses
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
    /// flag to signal that handle is set
    bool m_handleSet;
};
#endif