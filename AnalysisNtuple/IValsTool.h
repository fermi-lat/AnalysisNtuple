/** @file IValsTool.h
@brief common abstract inteface for all the XxxValsTools
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/AnalysisNtuple/IValsTool.h,v 1.19 2008/02/28 21:42:51 lsrea Exp $
*/

#ifndef _H_IValsTool
#define _H_IValsTool

#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/MsgStream.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IValsTool("IValsTool", 6 , 0); 

/** @class IValsTool
* @brief Abstract interface for the XxxValsTools, including visitor
*
* @author Leon Rochester
*
*/

namespace {
    enum check { NOCALC = -1, CALC = 0, CHECK = 1};
    }

class   IValsTool : virtual public IAlgTool
{
public:
    
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IValsTool; }
    
    /// get a particular value, using ntuple name
    virtual StatusCode getVal(std::string varName, double& value, int check = CALC) = 0;
    virtual StatusCode getVal(std::string varName, float& value, int check = CALC) = 0;
    virtual StatusCode getVal(std::string varName, int& value, int check = CALC) = 0;
    virtual StatusCode getVal(std::string varName, unsigned int& value, int check = CALC) = 0;
    /// get a particular value, using ntuple name, with calc checking
    virtual StatusCode getValCheck(std::string varName, double& value) =0;
    virtual StatusCode getValCheck(std::string varName, float& value) =0;
    virtual StatusCode getValCheck(std::string varName, int& value) =0;
    virtual StatusCode getValCheck(std::string varName, unsigned int& value) =0;
    /// output the names and values, either all (default) or just one;
    virtual StatusCode browse(MsgStream log, const std::string varName = "") =0;
    /// let the user trigger her own calculation
    virtual StatusCode doCalcIfNotDone() = 0;
    /// number of times a tool did its calculation for this event
    virtual int getCalcCount() = 0;
    /// return a string containing the value, as above
    virtual StatusCode getVal(std::string varName, std::string& value, int check = CALC) = 0;
    virtual StatusCode getValCheck(std::string varName, std::string& value) =0;
    virtual void setLoadOrder(int index) = 0;
    virtual int  getLoadOrder() = 0;
    virtual bool isLoaded() = 0;
    
    /** @class Visitor 
    @brief calls the client successively with the names and (ref to) values

    See NtupleVisitor for an example.
    */
    class Visitor 
    {
    public:
    /// visitor callback controls further actions of server
    enum eVisitorRet {
        /// normal return: continue processing 
        CONT,        
        /// client has all information desired; no more traversal 
        USER_DONE, 
        /// client has serious error; abort 
        ERROR,
        /// not used by client. Returned by traverse at end of normal processing
        DONE  
    };
        
        /// callback to send varnames and values to the client
        virtual Visitor::eVisitorRet analysisValue(std::string varName,
            const double& value) const =0;
        virtual Visitor::eVisitorRet analysisValue(std::string varName,
            const float& value) const =0;
        virtual Visitor::eVisitorRet analysisValue(std::string varName,
            const int& value) const =0;
        virtual Visitor::eVisitorRet analysisValue(std::string varName,
            const unsigned int& value) const =0;
        virtual Visitor::eVisitorRet analysisValue(std::string varName,
            const char* value) const =0;
    }; 
    
    
    /// sets up callback method for user to access the data
    virtual Visitor::eVisitorRet traverse(IValsTool::Visitor* v, 
        const bool checkCalc=true) = 0;   
};

#endif  // _H_IValsTool
