// $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/test/test_AnalysisNtuple.cxx,v 1.4 2003/02/26 01:29:07 lsrea Exp $

// Include files
// Gaudi system includes

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/AlgTool.h"


#include "AnalysisNtuple/IValsTool.h"

#include <vector>

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** 
A simple algorithm.
*/

class test_AnalysisNtuple : public Algorithm {
public:
    test_AnalysisNtuple(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    
    //! number of times called
    
    int m_count; 
    
    //tool stuff
    IToolSvc * m_pToolSvc;
    std::vector<IValsTool*> m_toolvec;
};

//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_AnalysisNtuple );

static const AlgFactory<test_AnalysisNtuple>  Factory;
const IAlgFactory& test_AnalysisNtupleFactory = Factory;

//------------------------------------------------------------------------
//! ctor

test_AnalysisNtuple::test_AnalysisNtuple(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
    
}

//------------------------------------------------------------------------

//! set parameters and attach to various perhaps useful services.

StatusCode test_AnalysisNtuple::initialize(){
    StatusCode sc   = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // set up tools
    m_pToolSvc = 0;
    
    sc = service("ToolSvc", m_pToolSvc, true);
    if (!sc.isSuccess ()){
        log << MSG::INFO << "Can't find ToolSvc, will quit now" << endreq;
        return StatusCode::FAILURE;
    }

    const char * toolnames[] = {"CalValsTool","TkrValsTool",
    "McValsTool","TkrHitValsTool", "VtxValsTool", "GltValsTool", "AcdValsTool"};

    for( int i =0; i< 7; ++i){
        m_toolvec.push_back(0);
        sc = m_pToolSvc->retrieveTool(toolnames[i], m_toolvec.back());
        if( sc.isFailure() ) {
            log << MSG::ERROR << "Unable to find a  tool" << toolnames[i]<< endreq;
            return sc;
        }
    }
  
    return sc;
}

//------------------------------------------------------------------------

//! process an event
StatusCode test_AnalysisNtuple::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    double answer;
    
    //do a browse
    m_toolvec[5]->browse();
    
    // check browse() against getVal() for each tool

    log << MSG::INFO << "Checks of getVal() against browse():" << endreq << endreq;
    
    m_toolvec[0]->browse("CAL_EneSum_Corr");
    sc = m_toolvec[0]->getVal("CAL_EneSum_Corr", answer);
    std::cout << "  compared to: " << answer << std::endl;
    
    m_toolvec[1]->browse("TKR_Energy_Corr");
    sc = m_toolvec[1]->getVal("TKR_Energy_Corr", answer);
    std::cout << "  compared to: " << answer << std::endl;
    
    m_toolvec[4]->browse("VTX_zdir");
    sc = m_toolvec[4]->getVal("VTX_zdir", answer);
    std::cout << "  compared to: " << answer << std::endl;
    
    m_toolvec[6]->browse("ACD_TileCount");
    sc = m_toolvec[6]->getVal("ACD_TileCount", answer);
    std::cout << "  compared to: " << answer << std::endl;
    
    m_toolvec[2]->browse("MC_x_err");
    sc = m_toolvec[2]->getVal("MC_x_err", answer);
    std::cout << "  compared to: " << answer << std::endl;
    
    m_toolvec[5]->browse("TRG_zDir");
    sc = m_toolvec[5]->getVal("TRG_zDir", answer);
    std::cout << "  compared to: " << answer << std::endl;
    
    m_toolvec[3]->browse("TKR_Hits_In_Lyr_12");
    sc = m_toolvec[3]->getVal("TKR_Hits_In_Lyr_12", answer);
    std::cout << "  compared to: " << answer << std::endl;
    
    
    return sc;
}



//------------------------------------------------------------------------

//! clean up, summarize
StatusCode test_AnalysisNtuple::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
    return sc;
}






