// $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/test/test_AnalysisNtuple.cxx,v 1.1.1.1 2003/02/23 21:31:10 burnett Exp $

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
    IValsTool* m_caltool;
    IValsTool* m_tkrtool;
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
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

        
    // set up tools
    m_pToolSvc = 0;
    
    sc = service("ToolSvc", m_pToolSvc, true);
    if (!sc.isSuccess ()){
        log << MSG::INFO << "Can't find ToolSvc, will quit now" << endreq;
        return StatusCode::FAILURE;
    }


    
    sc = m_pToolSvc->retrieveTool("CalValsTool", m_caltool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find a  tool" << endreq;
        return sc;
    }
    
    
    sc = m_pToolSvc->retrieveTool("TkrValsTool", m_tkrtool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find a  tool" << endreq;
        return sc;
    }

    return sc;
}

//------------------------------------------------------------------------

//! process an event
StatusCode test_AnalysisNtuple::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << " No test yet!" << endreq;

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







