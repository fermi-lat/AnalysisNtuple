
/** @file AcdValsTool.cxx
@brief Calculates the Adc analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/AcdValsTool.cxx,v 1.12 2004/08/23 19:56:44 heather Exp $
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
//#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Digi/AcdDigi.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
// Point used by AcdDigi
//#include "CLHEP/Geometry/Point3D.h"  <=== check this
// Point used by TKR
//#include "geometry/Point.h"  <=== check this

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
    
    void tkrHitsCount();

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
    double ACD_TkrHitsCountTop;
    double ACD_TkrHitsCountRows[4];

    IGlastDetSvc *m_detSvc;
    double m_vetoThresholdMeV;
    double m_tkrHitsCountCut;
    
    
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

    // in mm
    declareProperty("tkrHitsCountCut", m_tkrHitsCountCut=250.0);
}

StatusCode AcdValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

   if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;
   
   setProperties();
    // get the services
    
   // use the pointer from ValBase
   
   if( serviceLocator() ) {
                
        // find GlastDevSvc service
        if (service("GlastDetSvc", m_detSvc, true).isFailure()){
            log << MSG::INFO << "Couldn't find the GlastDetSvc!" << endreq;
            log << MSG::INFO << "Will be unable to calculate ACD_TkrHitsCount" << endreq;
            m_vetoThresholdMeV = 0.4;
        } else {
            StatusCode sc = m_detSvc->getNumericConstByName("acd.vetoThreshold", &m_vetoThresholdMeV);
            if (sc.isFailure()) {
                log << MSG::INFO << "Unable to retrieve threshold, setting the value to 0.4 MeV" << endreq;
                m_vetoThresholdMeV = 0.4;
   
             }
         }
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

    addItem("AcdTkrHitsCountTop", &ACD_TkrHitsCountTop);
    addItem("AcdTkrHitsCountR0", &ACD_TkrHitsCountRows[0]);
    addItem("AcdTkrHitsCountR1", &ACD_TkrHitsCountRows[1]);
    addItem("AcdTkrHitsCountR2", &ACD_TkrHitsCountRows[2]);
    addItem("AcdTkrHitsCountR3", &ACD_TkrHitsCountRows[3]);

    zeroVals();
    
    return sc;
}


StatusCode AcdValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    tkrHitsCount();
    
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

void AcdValsTool::tkrHitsCount() {

    // Purpose and Method:  Count the number of TkrClusters within some distance
    //  of hit ACD tiles.
    // Loops over all AcdDigis (skips those that are not "hit") 
    // Then for each TkrCluster, determine the distance between the center of 
    // The ACD tile and the TkrCluster.  If the distance is below m_tkrHitsCountCut
    // Then count this hit, either as top, R0, R1, R2, R3, depending upon what
    // type of ACD tile this is.

    MsgStream log(msgSvc(), name());
    if (!m_detSvc) return;

    SmartDataPtr<Event::AcdDigiCol> acdDigiCol(m_pEventSvc, EventModel::Digi::AcdDigiCol);
    if (!acdDigiCol) return;

    SmartDataPtr<Event::TkrClusterCol> pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    if (!pClusters) return;


    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = acdDigiCol->begin(); acdDigiIt != acdDigiCol->end(); acdDigiIt++) 
    {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        if (!acdId.tile()) continue;  // skip ribbons

        // Use MC Energy cut for now - until we move to veto discrim
        if ((*acdDigiIt)->getEnergy() < m_vetoThresholdMeV) continue; 
        // Try using the veto discriminator
        // If neither veto discrim from either PMT is high - we should skip
        //if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && 
        //     (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue;
        
        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
        std::string str;
        std::vector<double> dim;
        StatusCode sc = m_detSvc->getShapeByID(volId, &str, &dim);
        if ( sc.isFailure() ) {
            log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
            continue;
        }
        HepTransform3D transform;
        sc = m_detSvc->getTransform3DByID(volId, &transform);
        if (sc.isFailure() ) {
            log << MSG::WARNING << "Failed to get transformation" << endreq;
            continue;
        }
		
        HepPoint3D center(0., 0., 0.);
        HepPoint3D acdCenter = transform * center;

        // Loop over the clusters and calculate the distance from cluster
        // to ACD tile center
        int numClusters = pClusters->nHits();
        int icluster;
        for (icluster = 0; icluster < numClusters; icluster++) {
            Event::TkrCluster *clusterTds = pClusters->getHit(icluster);
            Point clusterPos = clusterTds->position();
            HepPoint3D clusterP(clusterPos.x(), clusterPos.y(), clusterPos.z());
            double dist = clusterP.distance(acdCenter);
            if (dist < m_tkrHitsCountCut) {
                if (acdId.top()) ++ACD_TkrHitsCountTop;
                if (acdId.side()) {
                    unsigned int row = acdId.row();
                    ++ACD_TkrHitsCountRows[row];
                }
            }
        }
    }


    return;
}

