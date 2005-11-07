/** @file McValsTool.cxx
@brief Calculates the Mc analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/McValsTool.cxx,v 1.24 2005/10/26 00:05:50 lsrea Exp $
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
#include "Event/TopLevel/MCEvent.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

// Reconstructed Tracks.... 
#include "Event/Recon/TkrRecon/TkrTrack.h"
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
    
    //Attempt to calculate the energy exiting the tracker
    double getEnergyExitingTkr(Event::McParticle* mcPart);
    
    //Pure MC Tuple Items
    float MC_SourceId;
    float MC_Id;
    float MC_Charge;
    float MC_Energy;
    float MC_LogEnergy;
    float MC_EFrac;
    float MC_OpenAngle; 
    float MC_TkrExitEne;

    
    float MC_x0;
    float MC_y0;
    float MC_z0;
    
    float MC_xdir;
    float MC_ydir;
    float MC_zdir;
    
    //MC - Compared to Recon Items
    float MC_x_err; 
    float MC_y_err;
    float MC_z_err;
    
    float MC_xdir_err; 
    float MC_ydir_err;
    float MC_zdir_err;
    
    float MC_dir_err;
    float MC_TKR1_dir_err;
    float MC_TKR2_dir_err;

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

    /** @page anatup_vars 
    @section mcvalstool McValsTool Variables

<table>
<tr><th> Variable </th><th> Description
<tr><td> McSourceId 
<td>        Unique integer associated with each MC source type; 
            from McEvent header replaces Mc_src_Id in merit ntuple 
<tr><td> McId 
<td>        StdHepId of primary (-13 = mu+, 22 = gamma, etc.) 
<tr><td> McCharge 
<td>        Charge of primary 
<tr><td> McEnergy 
<td>        Kinetic energy of the generated primary particle 
<tr><td> McLogEnergy 
<td>        log10(McEnergy) 
<tr><td> McEFrac 
<td>        Fraction of incident energy in highest-energy daughter 
<tr><td> McOpeningAngle 
<td>        Actual opening angle between the first and second daughters 
            of the promary as generated, (For a primary photon, 
            these will ordinarily be the electron and positron.) 
<tr><td> McTkrExitEne 
<td>        Attempt to calculate the total energy <strong>leaving</strong> the tracker volume 
<tr><td> Mc[X/Y/Z]0 
<td>        [x/y/z] coordinate of photon conversion or charged particle origin 
<tr><td> Mc[X/Y/Z]Dir 
<td>        [x/y/z] direction cosine of primary particle 
<tr><td> Mc[X/Y]Err 
<td>        [x/y] (found) - [x/y] (Mc) (Mc position taken at the z of the 
            found vertex or first hit)
<tr><td> McZErr 
<td>        z(<strong>actual</strong> vertex or first hit) - McZ0 
<tr><td> Mc[X/Y/Z]DirErr 
<td>        [x/y/z]dir (found) - [x/y/z]dir (Mc )
<tr><td> McDirErr 
<td>        Angle between found direction and Mc direction (radians )
<tr><td> McTkr[1/2]DirErr 
<td>        Angle between direction of [best/second] track and Mc direction (radians) 
</table>
    */


    addItem("McSourceId",     &MC_SourceId);
    addItem("McId",           &MC_Id);  
    addItem("McCharge",       &MC_Charge);
    addItem("McEnergy",       &MC_Energy);  
    addItem("McLogEnergy",    &MC_LogEnergy);
    addItem("McEFrac",        &MC_EFrac);
    addItem("McOpenAngle",    &MC_OpenAngle);
    addItem("McTkrExitEne",   &MC_TkrExitEne);
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
    SmartDataPtr<Event::TkrTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol); 
    SmartDataPtr<Event::TkrVertexCol>       pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    // Recover MC Pointer
    SmartDataPtr<Event::McParticleCol> pMcParticle(m_pEventSvc, EventModel::MC::McParticleCol);
    SmartDataPtr<Event::MCEvent> pMcEvent(m_pEventSvc, EventModel::MC::Event);
    
    if(pMcEvent) {
        MC_SourceId = pMcEvent->getSourceId();
    }
    
    MC_Energy = -1;
    
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

        //Attempt to estimate energy exiting the tracker
        MC_TkrExitEne = getEnergyExitingTkr(*pMCPrimary);

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
                Vector Mc_t1 = Vector(Mc_p1.x(),Mc_p1.y(), Mc_p1.z()).unit();
                Vector Mc_t2 = Vector(Mc_p2.x(),Mc_p2.y(), Mc_p2.z()).unit();
                double dot_prod = Mc_t1*Mc_t2;
                if(dot_prod > 1.) dot_prod = 1.;
                MC_OpenAngle = acos(dot_prod);
            }  
        }

        // This should really be a test on vertices, not tracks
        if(!pTracks) return sc; 
        int num_tracks = pTracks->size(); 
        if(num_tracks <= 0 ) return sc;
        
        // Get track energies and event energy
        Event::TkrTrackColConPtr pTrack1 = pTracks->begin();
        const Event::TkrTrack*   track_1 = *pTrack1;
        
        double e1 = track_1->getInitialEnergy();
        double gamEne = e1; 
        double e2 = 0.; 
        if(num_tracks > 2) {
            pTrack1++;
            const Event::TkrTrack* track_2 = *pTrack1;
            e2 = track_2->getInitialEnergy();
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
            
            SmartRefVector<Event::TkrTrack>::const_iterator pTrack1 = gamma->getTrackIterBegin();  
            const Event::TkrTrack* track_1 = *pTrack1;
            
            Point  x1 = track_1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
            Vector t1 = track_1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
            double cost1tMC = t1*Mc_t0;
            
            MC_TKR1_dir_err  = acos(cost1tMC);
            
            if(nParticles > 1) {
                pTrack1++;
                const Event::TkrTrack* track_2 = *pTrack1;
                Point  x2 = track_2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
                Vector t2 = track_2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
                double cost2tMC = t2*Mc_t0;
                
                MC_TKR2_dir_err  = acos(cost2tMC);
            }
        } 
    }
    return sc;
}

double McValsTool::getEnergyExitingTkr(Event::McParticle* mcPart)
{
    //This function attempts to estimate the energy which exist the tracker during 
    //an event. The idea is to recursively traverse the McParticle tree and look for
    //particles which start in the tracker (so, for example, the electron and positron
    //from a gamma conversion) and then exit, keep track of the energy of these particles.
    //Due to the vagaries of our simulation, several interesting problems can arise and
    //an attempt is made to deal with these (see notes below).
    //A flaw here is that NO attempt is made to account for continuous energy loss of 
    //charged particles in the tracker. Hopefully, a more sophisticated routine will 
    //appear soon.
    double tkrExitEne = 0.;
    double partEnergy = 0.;
    double partLostE  = 0.;

    //idents::VolumeIdentifier initial = mcPart->getInitialId();
    idents::VolumeIdentifier final   = mcPart->getFinalId();

    //This should check that the initial point of the particle is in the tracker
    //and that its final point is outside. So, it started in the tracker and 
    //carried its energy outside (could be anywhere)
    ////if (initial.size() == 9 && initial[0] == 0 && initial[3] == 1 &&
    ////    final.size() == 9 && (final[0] != 0 || final[3] != 1)       )
    if (final.size() == 9 && (final[0] != 0 || final[3] != 1) )
    {
        //Set total initial energy of this particle...
        partEnergy = mcPart->initialFourMomentum().t();
    }

    //Next step is to go through the daughter list to find the tracker exiting
    //energy due to them
    SmartRefVector<Event::McParticle> daughters = mcPart->daughterList();

    for(SmartRefVector<Event::McParticle>::iterator partIter = daughters.begin(); 
                                                    partIter < daughters.end(); 
                                                    partIter++)
    {
        Event::McParticle* daughter = *partIter;

        //First we check to see if the daughter was created in the tracker. 
        idents::VolumeIdentifier daughterInitial = daughter->getInitialId();

        if (daughterInitial.size() == 9 && daughterInitial[0] == 0 && daughterInitial[3] ==1)
        {
            //Keep track of the energy this daughter takes from its parent
            partLostE += daughter->initialFourMomentum().t();

            //Now get the energy exiting the tracker due to this particle
            tkrExitEne += getEnergyExitingTkr(daughter);
        }
    }

    //Ok, now correct the mcPart energy for that it lost traversing the tracker by 
    //creating daughter particles
    partEnergy -= partLostE;

    if (partEnergy < 0.) partEnergy = 0.;

    //Finally, add this to the energy exiting the tracker
    tkrExitEne += partEnergy;

    //This should now be the energy exiting the tracker due to those particles created 
    //in the tracker. 
    return tkrExitEne;
}
