#ifndef ValBase_h
#define ValBase_h

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"


#include "Event/TopLevel/Event.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "AnalysisNtuple/IValsTool.h"
#include <string>
#include <map>

// Some def's and functions to be used later
class ValBase : public IValsTool, virtual public IIncidentListener
{
public:

typedef std::map<std::string, double*> valMap;
typedef valMap::iterator mapIter;

    ValBase() : m_newEvent(true), m_handleSet(false) {m_ntupleMap.clear();}
    ~ValBase() {return;}

    /// clear map values
    virtual void zeroVals();
    /// fill ntuple values into tupleSvc
    virtual StatusCode fillNtuple(INTupleWriterSvc* pSvc, std::string tupleName);
    /// put one value into the ntuple
    virtual StatusCode addNtupleValue(std::string valName, INTupleWriterSvc* pSvc,
        std::string tupleName);
    /// check if calculation is already done for this event
    virtual StatusCode doCalcIfNotDone();
    /// get a particular value, using ntuple name
    virtual StatusCode getVal(std::string varName, double& value);
    /// output the list of names
    virtual void announceBadName(std::string varName);
    /// output the names and values, either all (default) or just one;
    virtual void browseValues(std::string varName = "");
    /// store the data provider
    //virtual void setEventSvc(IDataProviderSvc* svc);
    /// pass in the pointer to the incident service
    virtual void setIncSvc(IIncidentSvc* incsvc);
    /// this is called by the incident service at the beginning of an event
    virtual void handle(const Incident& inc);

    /// calculate all values, over-ridden by XxxValsTool
    virtual StatusCode calculate();

protected:

    /// map containing ntuple names, and pointers to the ntuple variables
    valMap m_ntupleMap;
    /// flag to signal new event
    bool m_newEvent;
    /// flag to signal that handle is set
    bool m_handleSet;
    /// pointer to incident service
    IIncidentSvc* m_incSvc;
};
#endif