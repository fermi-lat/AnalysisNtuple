
#ifndef _H_IValsTool
#define _H_IValsTool

#include "GaudiKernel/IAlgTool.h"

#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IValsTool("IValsTool", 2 , 3); 

/** @class ITkrValsTool
* @brief Abstract interface for tracker ntuple values 
*
* @author Leon Rochester
* $Header$
*
*/

/** Visitor callbacks can indicate whether traversal should continue
       or not.
     CONT        normal return: continue processing
     USER_DONE   client has all information desired; no more traversal
     ERROR       client has serious error; abort
     DONE        not used by client.  Will be returned by 
                 BadStrips::traverse   in case processing was normal.
*/

class ValsVisitor 
{
public:
    enum eVisitorRet {CONT, USER_DONE, ERROR, DONE};

    /// callback to send varnames and values to the client
    virtual ValsVisitor::eVisitorRet analysisValue(std::string VarName,
        double& value) const =0;
}; 

class   IValsTool : virtual public IAlgTool
{
public:
    
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IValsTool; }
    
    /// get a particular value, using ntuple name
    virtual StatusCode getVal(std::string varName, double& value) =0;
    /// output the names and values, either all (default) or just one;
    virtual StatusCode browse(std::string varName = "") =0;

    /// sets up callback method for user to access the data
    virtual ValsVisitor::eVisitorRet traverse(ValsVisitor* v) = 0;   
    /// handle from incident service, called at the beginning of each event
    virtual void handle(const Incident& inc) = 0;    
};

#endif  // _H_ITkrValsTool
