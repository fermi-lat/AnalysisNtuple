/** @file McKludgeValsTool.cxx
@brief Inserts sentinel values for several variables needed by ClassifyAlg
@author Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/McKludgeValsTool.cxx,v 1.1 2007/10/02 16:26:45 lsrea Exp $
*/
// Include files

#include "ValBase.h"

// for the moment, this code doesn't need any external data
// so I've commented out most of the includes

#include "GaudiKernel/MsgStream.h"
//#include "GaudiKernel/IDataProviderSvc.h"
//#include "GaudiKernel/SmartDataPtr.h"
//#include "GaudiKernel/SmartDataLocator.h"
//#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
//#include "GaudiKernel/IParticlePropertySvc.h"
//#include "GaudiKernel/ParticleProperty.h"


/*! @class McKludgeValsTool
@brief calculates bogus Monte Carlo values

@authors Bill Atwood, Leon Rochester
*/

class McKludgeValsTool : public ValBase
{
public:
    
    McKludgeValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~McKludgeValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
       
    //Pure MC Tuple Items
    float MC_SourceId;
    float MC_Tkr1DirErr;
    float MC_Energy;
    float MC_DirErr;
    float MC_DirErrN;
    float MC_DirErrN1;

    // to decode the particle charge
    //IParticlePropertySvc* m_ppsvc;    
};

// Static factory for instantiation of algtool objects
static ToolFactory<McKludgeValsTool> s_factory;
const IToolFactory& McKludgeValsToolFactory = s_factory;

// Standard Constructor
McKludgeValsTool::McKludgeValsTool(const std::string& type, 
                       const std::string& name, 
                       const IInterface* parent)
                       : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

/** @page anatup_vars 
@section McKludgeValsTool McKludgeValsTool Variables
These variables are supplied as a temporary fix so that
ClassifyAlg doesn't crash when it finds them missing

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> McSourceId 
<td>F<td>   Unique integer associated with each MC source type; 
            set to -1
<tr><td> McEnergy 
<td>F<td>   set to -1
<tr><td> McDirErr 
<td>F<td>   set to -1
<tr><td> McTkr1DirErr 
<td>F<td>   set to -1 
<tr><td> McDirErrN 
<td>F<td>   set to -1 
<tr><td> McDirErrN1 
<td>F<td>   set to -1 
</table>
*/


StatusCode McKludgeValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;
    
    /* 
    if( serviceLocator() ) {
        if( service("ParticlePropertySvc", m_ppsvc, true).isFailure() ) {
            log << MSG::ERROR << "Service [ParticlePropertySvc] not found" << endreq;
        }
    } else {
        return StatusCode::FAILURE;
    }
    */
    
    // load up the map

    addItem("McSourceId",     &MC_SourceId);
    addItem("McEnergy",       &MC_Energy);  
    addItem("McDirErr",       &MC_DirErr);      
    addItem("McTkr1DirErr",   &MC_Tkr1DirErr); 
    addItem("McDirErrN",      &MC_DirErrN);      
    addItem("McDirErrN1",     &MC_DirErrN1);      
    
    zeroVals();
    
    return sc;
}


StatusCode McKludgeValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MC_SourceId   = -1.0f;
    MC_Energy     = -1.0f;
    MC_DirErr     = -1.0f;
    MC_Tkr1DirErr = -1.0f;
    MC_DirErrN    = -1.0f;
    MC_DirErrN1   = -1.0f;
   
    return sc;
}

