
#ifndef _H_IValsTool
#define _H_IValsTool

#include "GaudiKernel/IAlgTool.h"


// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IValsTool("IValsTool", 2 , 0); 

class INTupleWriterSvc;
/** @class ITkrValsTool
* @brief Abstract interface for tracker ntuple values 
*
* @author Leon Rochester
* $Header$
*
*/

class   IValsTool : virtual public IAlgTool
{
public:
 
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IValsTool; }
        /// clear map values
    //virtual void zeroVals();
    /// fill ntuple values into tupleSvc
    virtual StatusCode fillNtuple(INTupleWriterSvc* pSvc, std::string tupleName)=0;
    /// put one value into the ntuple
  //  virtual StatusCode addNtupleValue(std::string valName, INTupleWriterSvc* pSvc,
//        std::string tupleName);
    /// check if calculation is already done for this event
    //virtual StatusCode doCalcIfNotDone();
    /// get a particular value, using ntuple name
    virtual StatusCode getVal(std::string varName, double& value)=0;
    /// output the list of names
    virtual void announceBadName(std::string varName)=0;
    /// output the names and values, either all (default) or just one;
    virtual void browseValues(std::string varName = "")=0;
    /// store the data provider
//?    virtual void setEventSvc(IDataProviderSvc* svc);

    /// calculate all values, over-ridden by XxxValsTool
    virtual StatusCode calculate()=0;

};

#endif  // _H_ITkrValsTool
