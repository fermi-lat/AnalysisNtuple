
#ifndef _H_IValsTool
#define _H_IValsTool

#include "GaudiKernel/IAlgTool.h"

#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IValsTool("IValsTool", 2 , 2); 

class INTupleWriterSvc;
/** @class ITkrValsTool
* @brief Abstract interface for tracker ntuple values 
*
* @author Leon Rochester
* $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/AnalysisNtuple/IValsTool.h,v 1.2 2003/02/24 17:55:45 lsrea Exp $
*
*/

class   IValsTool : virtual public IAlgTool
{
public:
    
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IValsTool; }
    
    /// fill ntuple values into tupleSvc
    virtual StatusCode fillNtuple(INTupleWriterSvc* pSvc, std::string tupleName)=0;
    /// put one value into the ntuple, for rolling your own ntuple!
    virtual StatusCode fillNtuple(INTupleWriterSvc* pSvc,
        std::string tupleName, std::string varName)=0;
    /// get a particular value, using ntuple name
    virtual StatusCode getVal(std::string varName, double& value)=0;
    /// output the list of names
    //virtual void announceBadName(std::string varName)=0;
    /// output the names and values, either all (default) or just one;
    virtual void browseValues(std::string varName = "")=0;
    
    /// calculate all values, over-ridden by XxxValsTool
    virtual StatusCode calculate()=0;
    
    /// handle from incident service
    virtual void handle(const Incident& inc) = 0;
    
};

#endif  // _H_ITkrValsTool
