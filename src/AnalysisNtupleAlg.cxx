
// $Header$

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

// visitor to do the ntuple (see IValsTool.h)
class NtupleVisitor : virtual public ValsVisitor
{
public:
    NtupleVisitor(INTupleWriterSvc* ntupleSvc=0, std::string ntupleName="") 
        : m_ntupleSvc(ntupleSvc), m_ntupleName(ntupleName) {}
    virtual ValsVisitor::eVisitorRet analysisValue(std::string varName, double& value) const;

private:
    INTupleWriterSvc* m_ntupleSvc;
    std::string m_ntupleName;
};

ValsVisitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName, double& value) const
{ 
    StatusCode sc;
    if (m_ntupleSvc) {
        sc = m_ntupleSvc->addItem(m_ntupleName.c_str(), varName.c_str(), value );
        if (sc.isFailure()) return ValsVisitor::ERROR;
    }    
    return ValsVisitor::CONT;
}

/** @class AnalysisNtuple
@brief fills the ntuple vfrom the XxxValsTools
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
    IValsTool* m_TkrVals;
    IValsTool* m_VtxVals;
    IValsTool* m_CalVals;
    IValsTool* m_AcdVals;
    IValsTool* m_McVals;
    IValsTool* m_GltVals;
    IValsTool* m_TkrHitVals;

    ValsVisitor* m_visitor;    
};

static const AlgFactory<AnalysisNtupleAlg>  Factory;
const IAlgFactory& AnalysisNtupleAlgFactory = Factory;

AnalysisNtupleAlg::AnalysisNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
    // declare properties with setProperties calls
    declareProperty("tupleName",  m_tupleName="");    
}

StatusCode AnalysisNtupleAlg::initialize(){
    StatusCode   sc = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    if( m_tupleName.empty()) {
        log << MSG::INFO << "tupleName property not set!  No ntuple output"<<endreq;
    }
    
    // calc tools
    
    sc = toolSvc()->retrieveTool("AcdValsTool", m_AcdVals);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve AcdValsTool" << endreq;
        return fail;
    }
    
    sc = toolSvc()->retrieveTool("TkrValsTool", m_TkrVals);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrValsTool" << endreq;
        return fail;
    }
    
    sc = toolSvc()->retrieveTool("VtxValsTool", m_VtxVals);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve VtxValsTool" << endreq;
        return fail;
    }
    
    sc = toolSvc()->retrieveTool("CalValsTool", m_CalVals);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve CalValsTool" << endreq;
        return fail;
    }
    
    sc = toolSvc()->retrieveTool("McValsTool", m_McVals);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve McValsTool" << endreq;
        return fail;
    }
    
    sc = toolSvc()->retrieveTool("GltValsTool", m_GltVals);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve GltValsTool" << endreq;
        return fail;
    }
    
    sc = toolSvc()->retrieveTool("TkrHitValsTool", m_TkrHitVals);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrHitValsTool" << endreq;
        return fail;
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
    
    return sc;
}

StatusCode AnalysisNtupleAlg::execute()
{
    StatusCode   sc = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;
    
    MsgStream   log( msgSvc(), name() );
    
    if (!m_tupleName.empty()) {

        if(m_ntupleSvc->addItem(m_tupleName.c_str(), "NumCalls", m_count).isFailure()) {
            log << MSG::ERROR << "AddItem failed" << endreq;
            return fail;
        }
        
        // fill the ntuples (triggers the calculations)
        if(m_McVals->traverse(m_visitor)==ValsVisitor::ERROR) {
            log << MSG::ERROR << "McVals traversal failed" << endreq;
            return fail;
        }
        if(m_GltVals->traverse(m_visitor)==ValsVisitor::ERROR) {
            log << MSG::ERROR << "AdcVals traversal failed" << endreq;
            return fail;
        }
        /*
        // not in Bill's ntuple output 
        if(m_TkrHitsVals->traverse(m_visitor)==ValsVisitor::ERROR) {
            log << MSG::ERROR << "TkrHitsVals traversal failed" << endreq;
            return fail;
        }
        */
        if(m_TkrVals->traverse(m_visitor)==ValsVisitor::ERROR) {
            log << MSG::ERROR << "TkrVals traversal failed" << endreq;
            return fail;
        }
        if(m_VtxVals->traverse(m_visitor)==ValsVisitor::ERROR) {
            log << MSG::ERROR << "VtxVals traversal failed" << endreq;
            return fail;
        }
        if(m_CalVals->traverse(m_visitor)==ValsVisitor::ERROR) {
            log << MSG::ERROR << "CalVals traversal failed" << endreq;
            return fail;
        }
        if(m_AcdVals->traverse(m_visitor)==ValsVisitor::ERROR) {
            log << MSG::ERROR << "AdcVals traversal failed" << endreq;
            return fail;
        }
 
        log << MSG::DEBUG;
        bool debugStuff = log.isActive();
        log << endreq;

        if (debugStuff) {
            double answer;

            //do a browse
            m_GltVals->browse();

            // check browse() against getVal() for each tool

            m_CalVals->browse("CAL_EneSum_Corr");
            sc = m_CalVals->getVal("CAL_EneSum_Corr", answer);
            log << MSG::DEBUG << "  compared to: " << answer << endreq;
            
            m_TkrVals->browse("TKR_Energy_Corr");
            sc = m_TkrVals->getVal("TKR_Energy_Corr", answer);
            log << MSG::DEBUG << "  compared to: " << answer << endreq;
            
            m_VtxVals->browse("VTX_DOCA");
            sc = m_VtxVals->getVal("VTX_DOCA", answer);
            log << MSG::DEBUG << "  compared to: " << answer << endreq;
            
            m_AcdVals->browse("ACD_TileCount");
            sc = m_AcdVals->getVal("ACD_TileCount", answer);
            log << MSG::DEBUG << "  compared to: " << answer << endreq;
            
            m_McVals->browse("MC_x_err");
            sc = m_McVals->getVal("MC_x_err", answer);
            log << MSG::DEBUG << "  compared to: " << answer << endreq;
            
            m_GltVals->browse("TRG_zDir");
            sc = m_GltVals->getVal("TRG_zDir", answer);
            log << MSG::DEBUG << "  compared to: " << answer << endreq;
            
            m_TkrHitVals->browse("TKR_Hits_In_Lyr_12");
            sc = m_TkrHitVals->getVal("TKR_Hits_In_Lyr_12", answer);
            log << MSG::DEBUG << "  compared to: " << answer << endreq;
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