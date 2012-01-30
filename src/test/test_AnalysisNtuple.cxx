/** @file test_AnalysisNtuple.cxx
@brief Placeholder: the actual test module is AnalysisNtupleAlg
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/test/test_AnalysisNtuple.cxx,v 1.9 2011/12/12 20:36:36 heather Exp $
*/

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
The test program for AnalysisNtuple is AnalysisNtupleAlg
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
};

//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_AnalysisNtuple );

//static const AlgFactory<test_AnalysisNtuple>  Factory;
//const IAlgFactory& test_AnalysisNtupleFactory = Factory;
DECLARE_ALGORITHM_FACTORY(test_AnalysisNtuple);

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
    return sc;
}

//------------------------------------------------------------------------

//! process an event
StatusCode test_AnalysisNtuple::execute()
{

    StatusCode  sc = StatusCode::SUCCESS;
    ++m_count;

    std::cout << "test code is in AnalysisNtupleAlg!" << std::endl;
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






