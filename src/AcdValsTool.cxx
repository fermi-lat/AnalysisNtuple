
/** @file AcdValsTool.cxx
@brief Calculates the Adc analysis variables
@author Bill Atwood, Leon Rochester

$Header$
*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"

/** @class AcdValsTool.cxx 
@brief Calculates Acd Values
*/

class AcdValsTool : public ValBase
{
public:
    
    AcdValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~AcdValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    
    //Global ACDTuple Items
    double ACD_Total_Energy;
    double ACD_Tile_Count; 
    double ACD_DOCA;
    double ACD_ActiveDist;
    double ACD_GammaDOCA; 
    
};

// Static factory for instantiation of algtool objects
static ToolFactory<AcdValsTool> s_factory;
const IToolFactory& AcdValsToolFactory = s_factory;

// Standard Constructor
AcdValsTool::AcdValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode AcdValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

   if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;
   
    // get the services
    
   // use the pointer from ValBase
   
   if( serviceLocator() ) {
                
   } else {
       return StatusCode::FAILURE;
   }
      
    // load up the map
    addItem("AcdTotalEnergy", &ACD_Total_Energy);
    addItem("AcdTileCount",    &ACD_Tile_Count);
    addItem("AcdDoca",         &ACD_DOCA);
    addItem("AcdActiveDist",   &ACD_ActiveDist);
    addItem("AcdGammaDoca",    &ACD_GammaDOCA);
    
    zeroVals();
    
    return sc;
}


StatusCode AcdValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    
    SmartDataPtr<Event::AcdRecon>           pACD(m_pEventSvc,EventModel::AcdRecon::Event);

    // Recover Track associated info. (not currently used 
    //SmartDataPtr<Event::TkrFitTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    //SmartDataPtr<Event::TkrVertexCol>     pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    // Recover pointer to ACD info  

    //Make sure we have valid ACD data
    if (pACD)
    {
        ACD_Total_Energy  = pACD->getEnergy();
        ACD_Tile_Count    = pACD->getTileCount(); 
        ACD_DOCA          = pACD->getDoca();
        ACD_ActiveDist    = pACD->getActiveDist();
        ACD_GammaDOCA     = pACD->getGammaDoca();
    } else {
        return StatusCode::FAILURE;
    }
    
    return sc;
}
