
/** @file AcdValsTool.cxx
@brief Calculates the Adc analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/AcdValsTool.cxx,v 1.10 2004/02/13 06:22:38 lsrea Exp $
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

#include <algorithm>
#include <numeric>

/** @class AcdValsTool
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
    double ACD_ActDistTop;
    double ACD_ActDistR0;
    double ACD_ActDistR1;
    double ACD_ActDistR2;
    double ACD_tileTopCount;
    double ACD_tileCount0;
    double ACD_tileCount1;
    double ACD_tileCount2;
    double ACD_tileCount3;
    double ACD_ribbon_ActiveDist;

    
};
namespace {

    // predicate to identify top, (row -1) or  side  (row 0-2)
    class acd_row { 
    public:
        acd_row(int row):m_row(row){}
        bool operator() ( std::pair<idents::AcdId ,double> entry){
            return m_row==-1? entry.first.face() == 0 : entry.first.row()==m_row;
        }
        int m_row;
    };
    // used by accumulate to get maximum for a given row
    class acd_max_energy { 
    public:
        acd_max_energy(int row):m_row(row), m_acd_row(row){}
        double operator()(double energy, std::pair<idents::AcdId ,double> entry) {
            return m_acd_row(entry)? std::max(energy, entry.second) : energy;
        }
        int m_row;
        acd_row m_acd_row;
    };
}
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

    addItem("AcdActDistTop",&ACD_ActDistTop);
    addItem("AcdActDistSideRow0",&ACD_ActDistR0);
    addItem("AcdActDistSideRow1",&ACD_ActDistR1);
    addItem("AcdActDistSideRow2",&ACD_ActDistR2);

    addItem("AcdNoTop",        &ACD_tileTopCount);
    addItem("AcdNoSideRow0",   &ACD_tileCount0);
    addItem("AcdNoSideRow1",   &ACD_tileCount1);
    addItem("AcdNoSideRow2",   &ACD_tileCount2);   
    addItem("AcdNoSideRow3",   &ACD_tileCount3);   
    addItem("AcdRibbonActDist", &ACD_ribbon_ActiveDist);

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
        ACD_ribbon_ActiveDist = pACD->getRibbonActiveDist();

        const std::vector<double> & adist = pACD->getRowActDistCol();
        ACD_ActDistTop = adist[0];
        ACD_ActDistR0 = adist[1];
        ACD_ActDistR1 = adist[2];
        ACD_ActDistR2 = adist[3];

        //Code from meritAlg.... 
        // get the map of energy vs tile id: have to construct from two parallel vectors
       float m_acd_tileCount[5];
       const std::vector<double> energies = pACD->getEnergyCol();
       const std::vector<idents::AcdId>& ids = pACD->getIdCol();
       std::vector<double>::const_iterator eit = energies.begin();

       std::map<idents::AcdId, double> emap;
       for( std::vector<idents::AcdId>::const_iterator idit = ids.begin(); 
           idit != ids.end() && eit !=energies.end(); ++idit, ++ eit){
           emap[*idit]=*eit;
       }
     
      // use acd_row predicate to count number of top tiles
      m_acd_tileCount[0] = std::count_if(emap.begin(), emap.end(), acd_row(-1) );

      // use acd_row predicate to count number of tiles per side row
      if(true)for( int row = 0; row<=3; ++row){ 
          m_acd_tileCount[row+1] = std::count_if(emap.begin(), emap.end(), acd_row(row) );
      }
      ACD_tileTopCount = m_acd_tileCount[0];
      ACD_tileCount0 = m_acd_tileCount[1];
      ACD_tileCount1 = m_acd_tileCount[2];
      ACD_tileCount2 = m_acd_tileCount[3];       
      ACD_tileCount3 = m_acd_tileCount[4];       

    } else {
        return StatusCode::FAILURE;
    }
    
    return sc;
}
