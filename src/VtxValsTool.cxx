
/** @file VtxValsTool.cxx
@brief Calculates the Vtx analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/VtxValsTool.cxx,v 1.13 2004/06/16 00:34:00 lsrea Exp $
*/

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

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "Doca.h"

namespace {
    double erf(double x) {
        double t = 1./(1.+.47047*x);
        double results = 1. -(.34802 - (.09587 -.74785*t)*t)*t*exp(-x*x);
        return results;
    }
    double thrshold(double x) {
        if(x < 0) return (.5*(1. + erf(-x)));
        else      return (.5*(1. - erf( x)));
    }
    
}

/*! @class VtxValsTool
@brief calculates Vtx values

@authors Bill Atwood, Leon Rochester
*/

class VtxValsTool : public ValBase
{
public:
    
    VtxValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~VtxValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    
    //Vertexing Items
    double VTX_xdir;
    double VTX_ydir;
    double VTX_zdir;
    double VTX_Phi;
    double VTX_Theta;
    double VTX_x0;
    double VTX_y0;
    double VTX_z0;
    double VTX_Angle;
    double VTX_DOCA;
    double VTX_Head_Sep;
    double VTX_S1;
    double VTX_S2;
    
    double VTX_DOCA_Wgt;  
    double VTX_HS_Wgt; 
    double VTX_S1_Wgt;     
    double VTX_T12_Wgt;
    double VTX_T2Q_Wgt;       
    double VTX_Total_Wgt;
    
};

// Static factory for instantiation of algtool objects
static ToolFactory<VtxValsTool> s_factory;
const IToolFactory& VtxValsToolFactory = s_factory;

// Standard Constructor
VtxValsTool::VtxValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode VtxValsTool::initialize()
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
    
    // Pair reconstruction
    addItem("VtxXDir",      &VTX_xdir);     
    addItem("VtxYDir",      &VTX_ydir);     
    addItem("VtxZDir",      &VTX_zdir);     
    addItem("VtxPhi",       &VTX_Phi);  
    addItem("VtxTheta",     &VTX_Theta);  
    addItem("VtxX0",        &VTX_x0);       
    addItem("VtxY0",        &VTX_y0);       
    addItem("VtxZ0",        &VTX_z0);       
    addItem("VtxAngle",     &VTX_Angle);    
    addItem("VtxDOCA",      &VTX_DOCA);     
    addItem("VtxHeadSep",   &VTX_Head_Sep); 
    
    addItem("VtxS1",        &VTX_S1);       
    addItem("VtxS2",        &VTX_S2);       
    
    addItem("VtxDocaWgt",  &VTX_DOCA_Wgt); 
    addItem("VtxHSWgt",    &VTX_HS_Wgt);   
    
    addItem("VtxS1Wgt",    &VTX_S1_Wgt);   
    addItem("VtxT12Wgt",   &VTX_T12_Wgt);  
    addItem("VtxT2QWgt",   &VTX_T2Q_Wgt);
    addItem("VtxTotalWgt", &VTX_Total_Wgt);
    
    
    zeroVals();
    
    return sc;
}


StatusCode VtxValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrTrackCol>  
        pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol);
    SmartDataPtr<Event::TkrVertexCol>     
        pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
   
    if(!pVerts) return sc; 
    
    // Recover total constrained energy
    if (!pTracks) return sc;
    int nParticles = pTracks->size();
    
    if(nParticles < 1) return sc;
    
    Event::TkrTrackColConPtr pTrack1 = pTracks->begin();
    const Event::TkrTrack* track_1 = *pTrack1;
    
    double gamEne = track_1->getInitialEnergy(); 
    
    if(nParticles > 2) {
        pTrack1++;
        const Event::TkrTrack* track_2 = *pTrack1;
        
        gamEne += track_2->getInitialEnergy();
    }
    
    //Make sure we have valid reconstructed data
    if (pVerts) {
        
        // Get the first Vertex - First track of first vertex = Best Track
        Event::TkrVertexConPtr pVtxr = pVerts->begin(); 
        Event::TkrVertex*   gamma = *pVtxr++; 
        SmartRefVector<Event::TkrTrack>::const_iterator pTrack1 = gamma->getTrackIterBegin(); 
        const Event::TkrTrack* track_1 = *pTrack1;
        
        nParticles = gamma->getNumTracks(); 
        Point  x0 = gamma->getPosition();
        Vector t0 = gamma->getDirection();
        VTX_xdir      = t0.x();
        VTX_ydir      = t0.y();
        VTX_zdir      = t0.z();

        // this code replaces atan and acos used before (also signs are consistent now!)
        VTX_Phi       = (-t0).phi();
        if (VTX_Phi<0.0) VTX_Phi += 2*M_PI;
        VTX_Theta     = (-t0).theta();

        VTX_x0        = x0.x();
        VTX_y0        = x0.y();
        VTX_z0        = x0.z();
        
        Point  x1 = track_1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
        Vector t1 = track_1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);

        double gamEne = track_1->getInitialEnergy();
        
        if(nParticles > 1) {
            pTrack1++;
            const Event::TkrTrack* track_2 = *pTrack1;
            
            Point  x2 = track_2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
            Vector t2 = track_2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
            
            // Put in a re-calc. of vertex quantities left out 
            Doca docaObj = Doca(Ray(x1,t1), Ray(x2,t2)); 
            
            VTX_Head_Sep = (x1-x2).magnitude(); 
            
            double cost1t2 = t1*t2; 
            double t1t2  = acos(cost1t2); 
            VTX_Angle = t1t2;
            VTX_DOCA  = docaObj.docaRay1Ray2();
            VTX_S1    = docaObj.arcLenRay1();
            VTX_S2    = docaObj.arcLenRay2();
            
            if(VTX_xdir == t1.x() && VTX_ydir==t1.y()) VTX_Angle = -.1; 
            
            double dist = VTX_DOCA;
            double s1   = VTX_S1;
            //double s2   = VTX_S2; 
            double root2 = sqrt(2.);
            double q2    = track_2->getQuality();
            const double minEne = 1.0;
            double thisEne = std::max(minEne, gamEne);
            VTX_DOCA_Wgt = thrshold((dist-50./thisEne)/root2);  
            VTX_S1_Wgt   = thrshold((s1-2.5)/(root2*1.5 * thisEne/100.));     
            VTX_T12_Wgt  = thrshold((.01*t1t2*thisEne -.2)/(root2*.2));
            VTX_T2Q_Wgt  = thrshold((40. - q2)/(root2*20.)); 
            VTX_HS_Wgt   = thrshold((VTX_Head_Sep - 2.)/(root2*2.));
            
            VTX_Total_Wgt= log(std::max(1.e-8, VTX_HS_Wgt*VTX_T12_Wgt*VTX_T2Q_Wgt));
        }
    }
    
    return sc;
}
