/** @file McValsTool.cxx
@brief Calculates the Mc analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/McValsTool.cxx,v 1.13 2003/10/17 21:03:11 lsrea Exp $
*/
// Include files


#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

// Reconstructed Tracks.... 
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

/*! @class McValsTool
@brief calculates Monte Carlo values

@authors Bill Atwood, Leon Rochester
*/

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
    double MC_Id;
    double MC_Charge;
    double MC_Energy;
	double MC_EFrac; 
    double MC_LogEnergy;
    
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

    // to decode the particle charge
    IParticlePropertySvc* m_ppsvc;    
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
        if( service("ParticlePropertySvc", m_ppsvc, true).isFailure() ) {
            log << MSG::ERROR << "Service [ParticlePropertySvc] not found" << endreq;
        }
    } else {
        return StatusCode::FAILURE;
    }
    
    // load up the map

    addItem("McId",           &MC_Id);  
    addItem("McCharge",       &MC_Charge);
    addItem("McEnergy",       &MC_Energy);  
    addItem("McLogEnergy",    &MC_LogEnergy);
	addItem("McEFrac",        &MC_EFrac);
    addItem("McX0",           &MC_x0);           
    addItem("McY0",           &MC_y0);           
    addItem("McZ0",           &MC_z0);           
    
    addItem("McXDir",         &MC_xdir);         
    addItem("McYDir",         &MC_ydir);         
    addItem("McZDir",         &MC_zdir);         
    
    addItem("McXErr",         &MC_x_err);        
    addItem("McYErr",         &MC_y_err);        
    addItem("McZErr",         &MC_z_err);        
    
    addItem("McXDirErr",      &MC_xdir_err);     
    addItem("McYDirErr",      &MC_ydir_err);     
    addItem("McZDirErr",      &MC_zdir_err);     
    
    addItem("McDirErr",       &MC_dir_err);      
    addItem("McTkr1DirErr",   &MC_TKR1_dir_err); 
    addItem("McTkr2DirErr",   &MC_TKR2_dir_err);   
    
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
    
    if (pMcParticle) {
        
        Event::McParticleCol::const_iterator pMCPrimary = pMcParticle->begin();
        // Skip the first particle... it's for bookkeeping.
        // The second particle is the first real propagating particle.
        pMCPrimary++;

        Event::McParticle::StdHepId hepid= (*pMCPrimary)->particleProperty();
        MC_Id = (double)hepid;
        ParticleProperty* ppty = m_ppsvc->findByStdHepID( hepid );
        if (ppty) {
            std::string name = ppty->particle(); 
            MC_Charge = ppty->charge();          
        }
        
        HepPoint3D Mc_x0;
        // launch point for charged particle; conversion point for neutral
        Mc_x0 = (MC_Charge==0 ? (*pMCPrimary)->finalPosition() : (*pMCPrimary)->initialPosition());
        HepLorentzVector Mc_p0 = (*pMCPrimary)->initialFourMomentum();
        // there's a method v.m(), but it does something tricky if m2<0
        double mass = sqrt(std::max(Mc_p0.m2(),0.0));
        
        Vector Mc_t0 = Vector(Mc_p0.x(),Mc_p0.y(), Mc_p0.z()).unit();
        
        //Pure MC Tuple Items
        MC_Energy = std::max(Mc_p0.t() - mass, 0.0);
        MC_LogEnergy = log10(MC_Energy);
        
        MC_x0     = Mc_x0.x();
        MC_y0     = Mc_x0.y();
        MC_z0     = Mc_x0.z();
        
        MC_xdir   = Mc_t0.x();
        MC_ydir   = Mc_t0.y();
        MC_zdir   = Mc_t0.z();

		if((*pMCPrimary)->daughterList().size() > 0) {
		    SmartRefVector<Event::McParticle> daughters = (*pMCPrimary)->daughterList();
			SmartRef<Event::McParticle> pp1 = daughters[0]; 
			std::string interaction = pp1->getProcess();
			if(interaction == "conv") { // Its a photon conversion; For comptons "compt" or brems "brem"  
				HepLorentzVector Mc_p1 = pp1->initialFourMomentum();
				SmartRef<Event::McParticle> pp2 = daughters[1];
				HepLorentzVector Mc_p2 = pp2->initialFourMomentum();
				double e1 = Mc_p1.t();
				double e2 = Mc_p2.t();
				MC_EFrac = e1/MC_Energy; 
				if(e1 < e2) MC_EFrac = e2/MC_Energy;
			}  
		}
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

			// Reference position errors at the start of recon track(s)
			double arc_len = (x0.z()-Mc_x0.z())/Mc_t0.z();
            HepPoint3D x_start = Mc_x0 + arc_len*Mc_t0;
            MC_x_err  = x0.x()-x_start.x(); 
            MC_y_err  = x0.y()-x_start.y();
            // except for z, use the difference between MC conversion point and start of vertex track
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