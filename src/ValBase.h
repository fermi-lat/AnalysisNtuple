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
typedef valMap::const_iterator constMapIter;

/// little class to preserve the order in the ntuple
class Order
{   
public:
    Order(std::string name="", mapIter iter=0) {}
    std::string name;
    mapIter iter;
    void setIter(mapIter newIter) { iter = newIter; }
    void setName(std::string newName) { name = newName; }
    std::string getName() { return name; }
    mapIter getIter() {return iter; }
};

    ValBase(const std::string& type, 
            const std::string& name, 
            const IInterface* parent);
    
    ~ValBase() 
    {
        for (int i=0; i<m_orderList.size(); i++) {
            Order* ord = m_orderList[i];
            delete ord;
        }
    }
    
    /// clear map values
    virtual void zeroVals();
    /// add an item to the map
    virtual void addItem(std::string varName, double* pValue);
    /// check if calculation is already done for this event
    virtual StatusCode doCalcIfNotDone();
    /// get a particular value, using ntuple name
    virtual StatusCode getVal(std::string varName, double& value);
    /// output the list of names
    virtual void announceBadName(std::string varName);
    /// output the names and values, either all (default) or just one;
    virtual StatusCode browse(std::string varName = "");
    /// this is called by the incident service at the beginning of an event
    virtual void handle(const Incident& inc);
    /// callback for visitor
    virtual ValsVisitor::eVisitorRet traverse(ValsVisitor* v);
    
    /// calculate all values, over-ridden by XxxValsTool
    virtual StatusCode calculate();

    // common initialization for subclasses
    virtual StatusCode initialize();

protected:
    
    /// map containing ntuple names, and pointers to the ntuple variables
    valMap m_ntupleMap;
    /// to preserve the order of the ntuple 
    //  probably should be a vector of pairs:  std::vector< std::pair<std::string, mapIter>* >;
    std::vector<Order*> m_orderList;
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