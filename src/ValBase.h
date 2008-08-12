/** @file ValBase.h
@brief header file for ValBase.cxx
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ValBase.h,v 1.37 2008/07/14 23:43:59 lsrea Exp $
*/

#ifndef ValBase_h
#define ValBase_h

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IIncidentListener.h"

#include "AnalysisNtuple/IValsTool.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>  // for M_PI, among others

class IIncidentSvc;
class IDataProviderSvc;

/** @class ValBase
@brief base class for the XxxValsTools

@author Leon Rochester

*/
namespace {

// LSR 14-Jul-08 code for ntuple types

    enum valType {DOUBLE, FLOAT, INT, UINT, ULONG64, STRING};
    class TypedPointer 
    {
    public:  
        TypedPointer(valType type, void* pointer, int dim=1) : m_type(type), 
            m_pointer(pointer), m_dim(dim)
        {}
        ~TypedPointer() {}

        valType getType()       { return    m_type; }
        void* getPointer()      { return    m_pointer; }
        int getDim()            { return    m_dim; }

// LSR 14-Jul-08 code for ntuple types

        void setVal(unsigned int val)   { 
            *(reinterpret_cast<unsigned int*>(getPointer())) = val; 
        }
        void setVal(unsigned long long val)   { 
            *(reinterpret_cast<unsigned long long*>(getPointer())) = val; 
        }
        void setVal(int val)    { *(reinterpret_cast<int*>(getPointer())) = val; }
        void setVal(float val)  { *(reinterpret_cast<float*>(getPointer())) = val; }
        void setVal(double val) { *(reinterpret_cast<double*>(getPointer())) = val; }

    private:
        valType m_type;
        void*   m_pointer;
        int     m_dim;
    };
}

class ValBase : public IValsTool,  public AlgTool,  virtual public IIncidentListener
{
public:
    typedef std::pair<std::string, TypedPointer*> valPair;
    typedef std::vector<valPair*> valMap;
    typedef valMap::iterator mapIter;
    typedef valMap::const_iterator constMapIter;

    ValBase(const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    ~ValBase(); 
    /// clear map values
    virtual void zeroVals();

// LSR 14-Jul-08 code for ntuple types
    /// add an item to the map
    virtual void addItem(std::string varName, double* pValue);
    virtual void addItem(std::string varName, float* pValue);
    virtual void addItem(std::string varName, int* pValue);
    virtual void addItem(std::string varName, unsigned int* pValue);
    virtual void addItem(std::string varName, unsigned long long* pValue);
    virtual void addItem(std::string varName, char* pValue);
    virtual StatusCode doCalcIfNotDone();
 
// LSR 14-Jul-08 code for ntuple types
   /// get a particular value, using ntuple name default forces calculation
    virtual StatusCode getVal(std::string varName, double& value, int check = CALC);
    virtual StatusCode getVal(std::string varName, float& value, int check = CALC);
    virtual StatusCode getVal(std::string varName, int& value, int check = CALC);
    virtual StatusCode getVal(std::string varName, unsigned int& value, int check = CALC);
    virtual StatusCode getVal(std::string varName, unsigned long long& value, int check = CALC);
    virtual StatusCode getVal(std::string varName, std::string& value, int check = CALC);
    /// get a particular value, using ntuple name, with calc checking (called by AnaTup)
 
// LSR 14-Jul-08 code for ntuple types
    virtual StatusCode getValCheck(std::string varName, double& value);
    virtual StatusCode getValCheck(std::string varName, float& value);
    virtual StatusCode getValCheck(std::string varName, int& value);
    virtual StatusCode getValCheck(std::string varName, unsigned int& value);
    virtual StatusCode getValCheck(std::string varName, unsigned long long& value);
    virtual StatusCode getValCheck(std::string varName, std::string& value);

    virtual bool getArrayArg(std::string varName, std::string& baseName,
        int& dim);
    virtual std::string getFullName(std::string varName, int dim);
   
    /// output the list of names
    virtual void announceBadName(std::string varName);
    /// output the names and values, either all (default) or just one;
    virtual StatusCode browse(MsgStream log, const std::string varName = "");
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

    /// AnaTup loaded this object
    virtual void setLoadOrder(int index) { m_loadOrder = index; }
    virtual bool isLoaded()              { return m_loadOrder>-1; }
    virtual int  getLoadOrder()          { return m_loadOrder; }
    
protected:
    StatusCode getTypedPointer(std::string varName, TypedPointer*& ptr, int check);

    /// some static methods

    /// sign of a number
    static double sign(double x) { return x>0 ? 1.: -1. ;}
    /// turn a global coordinate (tower, ladder, wafer) roughly into a local one
    static double globalToLocal(double x, double pitch, int n) {
        double xNorm = x/pitch + 0.5*n;
        return sign(x)*(fmod(fabs(xNorm),1.0) - 0.5)*pitch ;
    }

    static double circleFraction(double r) {
        double rl = (fabs(r) < 1.) ? fabs(r):1.; 
        double a_slice = 2.*(M_PI/4. - rl*sqrt(std::max(0.0,1.-rl*rl))/2. - asin(rl)/2.);
        double in_frac = 1.-a_slice/M_PI;
        if(r < 0.) in_frac = a_slice/M_PI;
        return in_frac;
    }

    static double circleFractionSimpson(double r, double angle_factor) {
        double slice_0 = circleFraction(r);
        double slice_p = circleFraction(r+angle_factor);
        double slice_m = circleFraction(r-angle_factor);
        return (slice_p + 4.*slice_0 + slice_m)/6.;
    }

    void printHeader(MsgStream& log);
    
    /// map containing ntuple names, and pointers to the ntuple variables
    valMap m_ntupleMap;
    /// pointer to incident service
    IIncidentSvc* m_incSvc;
    /// let ValBase handle the pointer to the data service, everyone uses it
    IDataProviderSvc* m_pEventSvc;
    /// flag to signal new event
    bool m_newEvent;
    /// flag to allow an always-calculate call if 0; if 1 checks and sets m_newEvent
    /// if -1 skips calculation
    int m_check;

    /// count calls to tools
    int m_calcCount;

    /// tells if this routine has been "loaded" by AnalysisNtupleAlg
    int m_loadOrder;

    /// Obvious "bad" value if an exception occurs whild computing output
    /// variables
    static const int s_badVal;


};
#endif
