// $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/src/test/test_AnalysisNtuple.cxx,v 1.2 2003/01/10 19:35:42 lsrea Exp $

// Include files
// Gaudi system includes

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "TkrUtil/ITkrFailureModeSvc.h"

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

  /// pointer to failure mode service
  ITkrFailureModeSvc* m_FailSvc;

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

    sc = service("TkrFailureModeSvc", m_FailSvc);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to find TkrFailureMode service" << endreq;
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







