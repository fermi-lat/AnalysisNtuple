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
class ValBase : public IValsTool {
public:

typedef std::map<std::string, double*> valMap;
typedef valMap::iterator mapIter;

    ValBase() : m_run(-1), m_event(-1) {m_ntupleMap.clear();}
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
    virtual void setEventSvc(IDataProviderSvc* svc);

    /// calculate all values, over-ridden by XxxValsTool
    virtual StatusCode calculate();


protected:

    /// map containing ntuple names, and pointers to the ntuple variables
    valMap m_ntupleMap;
    /// run number, to test first time through
    int m_run;
    /// event number, to test first time through
    int m_event;
    /// passed in by Tool
    IDataProviderSvc* pEventSvc;
};
#endif