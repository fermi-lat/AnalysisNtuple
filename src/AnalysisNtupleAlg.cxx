/** @file AnalysisNtupleAlg.cxx
@brief Uses the XxxValsTools to produce a comprehensive ntuple
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/AnalysisNtupleAlg.cxx,v 1.14 2003/04/01 19:44:54 lsrea Exp $
*/

// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h" 
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ISvcLocator.h"

// ntuple interface
#include "ntupleWriterSvc/INTupleWriterSvc.h"

// analysis tools
#include "AnalysisNtuple/IValsTool.h"

// for access to geometry perhaps
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

/** @class NtupleVisitor
@brief fields the callback from the tools; fills ntuple with names and values
@author Leon Rochester
*/
class NtupleVisitor : virtual public IValsTool::Visitor
{
public:
    /// constructor can set ntuple service and ntupleName
    NtupleVisitor(INTupleWriterSvc* ntupleSvc=0, std::string ntupleName="") 
        : m_ntupleSvc(ntupleSvc), m_ntupleName(ntupleName) {}
    virtual IValsTool::Visitor::eVisitorRet 
        analysisValue(std::string varName, const double& value) const;
    
private:
    /// pointer to the ntuple servic
    INTupleWriterSvc* m_ntupleSvc;
    /// name of the ntuple; should be the same as is set in NtupleWriterSvc
    std::string m_ntupleName;
};

IValsTool::Visitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName, 
                                                             const double& value) const
{ 
    StatusCode sc;
    if (m_ntupleSvc) {
        sc = m_ntupleSvc->addItem(m_ntupleName.c_str(), varName.c_str(), value );
        if (sc.isFailure()) return IValsTool::Visitor::ERROR;
    }    
    return IValsTool::Visitor::CONT;
}

/** @class AnalysisNtupleAlg
@brief fills the ntuple from the XxxValsTools
@author Leon Rochester
*/

class AnalysisNtupleAlg : public Algorithm {
public:
    AnalysisNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator);
    /// set parameters and attach to various perhaps useful services.
    StatusCode initialize();
    /// process one event
    StatusCode execute();
    /// clean up
    StatusCode finalize();
    
private: 
    /// number of times called
    int m_count; 
    
    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc *m_ntupleSvc;
    /// parameter to store the logical name of the ROOT file to write to
    std::string m_tupleName;
    
    /// Common interface to analysis tools
    std::vector<IValsTool*> m_toolvec;
    /// tool names
    std::vector<std::string> m_toolnames;
    
    IValsTool::Visitor* m_visitor;    
};

static const AlgFactory<AnalysisNtupleAlg>  Factory;
const IAlgFactory& AnalysisNtupleAlgFactory = Factory;

AnalysisNtupleAlg::AnalysisNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
    // declare properties with setProperties calls
    declareProperty("tupleName",  m_tupleName=""); 
    // so it looks like NTupleWriterSvc property, no harm having both!
    declareProperty("tuple_name",  m_tupleName="");   
    // List of tools to use -- maybe a bit kludgy, since the spelling need to be correct!
    declareProperty("toolList", m_toolnames);
    
}

StatusCode AnalysisNtupleAlg::initialize(){
    StatusCode   sc = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // calc tools - default is the full set!

    // Use the Job options service to set the Algorithm's parameters
    m_toolnames.clear();

    setProperties();

    //probably a better way to do this!
    // default set:
    std::string toolnames [] = {"Mc", "Glt", "TkrHit", "Tkr", "Vtx",  "Cal", "Acd", ""};

    int i;
    int namesSize;

    if (m_toolnames.empty()) {
        for (i=0; ; ++i) {
            if (toolnames[i]=="") break;
            m_toolnames.push_back(toolnames[i]+"ValsTool");
        }
    } else {
        namesSize = m_toolnames.size();
        for (i=0; i<namesSize; ++i) {
            m_toolnames[i] = m_toolnames[i]+"ValsTool";
        }
    }

    log << MSG::INFO << "Tools requested: ";
    namesSize = m_toolnames.size();
    for (i=0; i<namesSize; ++i) {
        log << m_toolnames[i] << " " ;
    }
    log << endreq;
    
    if( m_tupleName.empty()) {
        log << MSG::INFO << "tupleName property not set!  No ntuple output"<<endreq;
    }
    // set up tools
    IToolSvc* pToolSvc = 0;
    
    sc = service("ToolSvc", pToolSvc, true);
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Can't find ToolSvc, will quit now" << endreq;
        return StatusCode::FAILURE;
    }
        
    namesSize = m_toolnames.size();
    for (i =0; i<namesSize; ++i){
        m_toolvec.push_back(0);
        sc = pToolSvc->retrieveTool(m_toolnames[i], m_toolvec.back());
        if( sc.isFailure() ) {
            log << MSG::ERROR << "Unable to find tool: " << m_toolnames[i] << endreq;
            return sc;
        }
    }
        
    // get a pointer to our ntupleWriterSvc
    if (!m_tupleName.empty()) {
        if (service("ntupleWriterSvc", m_ntupleSvc, true).isFailure()) {
            log << MSG::ERROR 
                << "AnalysisNtupleAlg failed to get the ntupleWriterSvc" 
                << endreq;
            return fail;
        }
    }
    
    m_visitor = new NtupleVisitor(m_ntupleSvc, m_tupleName);
    
    log << MSG::DEBUG;
    if(log.isActive()) {
        log << "AnalysisNtuple called" << endreq;
    }
    log << endreq;

    if (!m_tupleName.empty()) {
                
        if(m_ntupleSvc->addItem(m_tupleName.c_str(), "NumCalls", 0).isFailure()) {
            log << MSG::ERROR << "AddItem failed" << endreq;
            return fail;
        }

        int size = m_toolvec.size();
        for( int i =0; i<size; ++i){
            if(m_toolvec[i]->traverse(m_visitor, false)==IValsTool::Visitor::ERROR) {
                log << MSG::ERROR << m_toolvec[i] << " traversal failed" << endreq;
                return fail;
            }
        }
    }

    log << MSG::INFO << "AnalysisNtuple set with " << m_toolvec.size() << " tools" << endreq;
    
    return sc;
}

StatusCode AnalysisNtupleAlg::execute()
{
    StatusCode   sc = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;
    
    MsgStream   log( msgSvc(), name() );

    /* test for missing first event
    if (m_count==0) {
        m_ntupleSvc->storeRowFlag(false);
        m_count++;
        return sc;
    }
    */
    
    if (!m_tupleName.empty()) {
        
        
        if(m_ntupleSvc->addItem(m_tupleName.c_str(), "NumCalls", m_count).isFailure()) {
            log << MSG::ERROR << "AddItem failed" << endreq;
            return fail;
        }
        ++m_count;
        
        int size = m_toolvec.size();
        for( int i =0; i<size; ++i){
            if(m_toolvec[i]->traverse(m_visitor)==IValsTool::Visitor::ERROR) {
                log << MSG::ERROR << m_toolvec[i] << " traversal failed" << endreq;
                return fail;
            }
        }
    }
    
    bool debugStuff = false;
    log << MSG::DEBUG;
    if (log.isActive()) {
        debugStuff = true;;
        log << "Debug display: ";
    }
    log << endreq;
    
    if(debugStuff) {
        double answer;
        
        //do a browse
        m_toolvec[0]->browse();

        int i;
        std::string varname;
        std::vector<std::string> varnames;
        varnames.clear();
        int namesSize = m_toolnames.size();
        for (i=0; i<namesSize; ++i) {
            std::string toolname = m_toolnames[i];
            if      (toolname=="McValsTool"     ) {varname = "McXErr";}
            else if (toolname=="GltValsTool"    ) {varname = "GltTotal";}
            else if (toolname=="TkrHitValsTool" ) {varname = "TkrHitsInLyr00";}
            else if (toolname=="TkrValsTool"    ) {varname = "TkrEnergyCorr";}
            else if (toolname=="VtxValsTool"    ) {varname = "VtxZDir";}
            else if (toolname=="CalValsTool"    ) {varname = "CalEneSumCorr";}
            else if (toolname=="AcdValsTool"    ) {varname = "AcdTileCount";}
            else                                  {varname = "";}
            varnames.push_back(varname);
        }
               
        // check browse() against getVal() for each tool
                
        int vecSize = m_toolvec.size();
        for(i=0; i<vecSize; ++i){
            varname = varnames[i];
            if (varname=="") continue;
            m_toolvec[i]->browse(varnames[i]);
            sc = m_toolvec[i]->getVal(varname, answer);
            log << MSG::DEBUG;
            if (log.isActive()) {
                log << "  compared to: " << answer;
            }
            log << endreq;    
        }        
    }     
    return sc;
}

StatusCode AnalysisNtupleAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
    return sc;
}