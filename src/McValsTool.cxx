// $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/MCValsTool.cxx,v 1.2 2003/02/26 01:12:43 lsrea Exp $

// Include files


#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

// Reconstructed Tracks.... 
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

class McValsTool : public ValBase
{
public:
    
    McValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~McValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    
    //Pure MC Tuple Items
    double MC_Energy;
    
    double MC_x0;
    double MC_y0;
    double MC_z0;
    
    double MC_xdir;
    double MC_ydir;
    double MC_zdir;
    
    //MC - Compared to Recon Items
    double MC_x_err; 
    double MC_y_err;
    double MC_z_err;
    
    double MC_xdir_err; 
    double MC_ydir_err;
    double MC_zdir_err;
    
    double MC_dir_err;
    double MC_TKR1_dir_err;
    double MC_TKR2_dir_err;
    
};

// Static factory for instantiation of algtool objects
static ToolFactory<McValsTool> s_factory;
const IToolFactory& McValsToolFactory = s_factory;

// Standard Constructor
McValsTool::McValsTool(const std::string& type, 
                       const std::string& name, 
                       const IInterface* parent)
                       : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode McValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;
  
    // get the services    
    
    if( serviceLocator() ) {               
    } else {
        return StatusCode::FAILURE;
    }
    
    // load up the map

    m_ntupleMap["MC_Energy"]        = &MC_Energy;       
    m_ntupleMap["MC_x0"]            = &MC_x0;           
    m_ntupleMap["MC_y0"]            = &MC_y0;           
    m_ntupleMap["MC_z0"]            = &MC_z0;           
    
    m_ntupleMap["MC_xdir"]          = &MC_xdir;         
    m_ntupleMap["MC_ydir"]          = &MC_ydir;         
    m_ntupleMap["MC_zdir"]          = &MC_zdir;         
    
    m_ntupleMap["MC_x_err"]         = &MC_x_err;        
    m_ntupleMap["MC_y_err"]         = &MC_y_err;        
    m_ntupleMap["MC_z_err"]         = &MC_z_err;        
    
    m_ntupleMap["MC_xdir_err"]      = &MC_xdir_err;     
    m_ntupleMap["MC_ydir_err"]      = &MC_ydir_err;     
    m_ntupleMap["MC_zdir_err"]      = &MC_zdir_err;     
    
    m_ntupleMap["MC_dir_err"]       = &MC_dir_err;      
    m_ntupleMap["MC_TKR1_dir_err"]  = &MC_TKR1_dir_err; 
    m_ntupleMap["MC_TKR2_dir_err"]  = &MC_TKR2_dir_err;   
    
    zeroVals();
    
    return sc;
}


StatusCode McValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>    pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    SmartDataPtr<Event::TkrVertexCol>       pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    // Recover MC Pointer
    SmartDataPtr<Event::McParticleCol> pMcParticle(m_pEventSvc, EventModel::MC::McParticleCol);
    
    //Make sure we have valid ACD data
    //Make sure we have valid cluster data
    if (pMcParticle) {
        
        // Get the first Track - it should be the "Best Track"
        Event::McParticleCol::const_iterator pMCTrack1 = pMcParticle->begin();
        
        HepPoint3D Mc_x0       = (*pMCTrack1)->finalPosition();
        HepLorentzVector Mc_p0 = (*pMCTrack1)->initialFourMomentum();
        
        Vector Mc_t0 = Vector(Mc_p0.x(),Mc_p0.y(), Mc_p0.z()).unit();
        
        //Pure MC Tuple Items
        MC_Energy = Mc_p0.t();
        
        MC_x0     = Mc_x0.x();
        MC_y0     = Mc_x0.y();
        MC_z0     = Mc_x0.z();
        
        MC_xdir   = Mc_t0.x();
        MC_ydir   = Mc_t0.y();
        MC_zdir   = Mc_t0.z();
        
        if(!pTracks) return sc; 
        int num_tracks = pTracks->size(); 
        if(num_tracks <= 0 ) return sc;
        
        // Get track energies and event energy
        Event::TkrFitConPtr pTrack1 = pTracks->begin();
        const Event::TkrFitTrackBase* trackBase = *pTrack1;
        const Event::TkrKalFitTrack* track_1 = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase);
        
        double e1 = track_1->getEnergy();
        double gamEne = e1; 
        double e2 = 0.; 
        if(num_tracks > 2) {
            pTrack1++;
            trackBase = *pTrack1;
            const Event::TkrKalFitTrack* track_2 = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase);
            e2 = track_2->getEnergy();
            gamEne += e2;
        }
        
        //Make sure we have valid reconstructed data
        if (pVerts) {
            
            // Get the first Vertex - First track of first vertex = Best Track
            Event::TkrVertexConPtr pVtxr = pVerts->begin(); 
            Event::TkrVertex*   gamma = *pVtxr++; 
            Point  x0 = gamma->getPosition();
            Vector t0 = gamma->getDirection();
            
            MC_x_err  = x0.x()-Mc_x0.x(); 
            MC_y_err  = x0.y()-Mc_x0.y();
            MC_z_err  = x0.z()-Mc_x0.z();
            
            MC_xdir_err = t0.x()-Mc_t0.x(); 
            MC_ydir_err = t0.y()-Mc_t0.y();
            MC_zdir_err = t0.z()-Mc_t0.z();
            
            double cost0tMC = t0*Mc_t0;
            
            MC_dir_err  = acos(cost0tMC);
            
            int nParticles = gamma->getNumTracks(); 
            
            SmartRefVector<Event::TkrFitTrackBase>::const_iterator pTrack1 = gamma->getTrackIterBegin();  
            trackBase = *pTrack1;
            SmartRef<Event::TkrKalFitTrack> track_1  = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase); 
            
            Point  x1 = track_1->getPosition();
            Vector t1 = track_1->getDirection();
            double cost1tMC = t1*Mc_t0;
            
            MC_TKR1_dir_err  = acos(cost1tMC);
            
            if(nParticles > 1) {
                pTrack1++;
                trackBase = *pTrack1;
                SmartRef<Event::TkrKalFitTrack> track_2 = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase);
                Point  x2 = track_2->getPosition();
                Vector t2 = track_2->getDirection();
                double cost2tMC = t2*Mc_t0;
                
                MC_TKR2_dir_err  = acos(cost2tMC);
            }
        } 
    }
    
    return sc;
}