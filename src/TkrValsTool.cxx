// $Header$

// Include files


#include "AnalysisNtuple/IValsTool.h"
#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GaudiKernel/IToolSvc.h"

//#include "TMAth.h"
 
#ifndef M_PI
#define M_PI = 3.14159265358979323846
#endif

double sign(double x) { return (x<0.) ? -1.:1.;} 

double twrEdgeT(double x, double y, double pitch, int &XY, int &outter) {
    double edge = 0.; 
    double x_twr = sign(x)*(fmod(fabs(x),pitch) - pitch/2.);
    double y_twr = sign(y)*(fmod(fabs(y),pitch) - pitch/2.);
    
    outter = 0; 
    
    if(fabs(x_twr) > fabs(y_twr)) {
        edge = pitch/2. - fabs(x_twr);
        XY = 1; 
        if(fabs(x) > 1.5*pitch) outter = 1;
    }
    else {
        edge = pitch/2. - fabs(y_twr);
        XY = 2;
        if(fabs(y) > 1.5*pitch) outter = 1;
    }
    return edge;
}

double circle_fracT(double r) {
    double rl = (fabs(r) < 1.) ? fabs(r):1.; 
    double a_slice = 2.*(M_PI/4. - rl*sqrt(1.-rl*rl)/2. - asin(rl)/2.);
    double in_frac = 1.-a_slice/M_PI;
    if(r < 0.) in_frac = a_slice/M_PI;
    return in_frac;
}    


class TkrValsTool : public AlgTool, public ValBase {
public:
    
    TkrValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~TkrValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    
    // some pointers to services
    
    /// pointer to tracker geometry
    ITkrGeometrySvc*       pTkrGeoSvc;
    /// pointer to event data service
    IDataProviderSvc*      m_pEventSvc;
    /// pointer to queryclusterstool
    ITkrQueryClustersTool* pQueryClusters;
    
    Event::TkrClusterCol*  m_Clusters; 
    Event::TkrFitTrackCol* m_Tracks;
    
    IKalmanParticle*       pKalParticle;
    
    
    //Global Track Tuple Items
    double Tkr_No_Tracks;
    double Tkr_Sum_KalEne; 
    double Tkr_Sum_ConEne;
    double Tkr_Energy;
    double Tkr_Energy_Sum;
    double Tkr_Energy_Corr;
    double Tkr_Edge_Corr;
    double Tkr_Total_Hits;
    double Tkr_Thin_Hits;
    double Tkr_Thick_Hits;
    double Tkr_Blank_Hits;  
    double Tkr_RadLength; 
    double Tkr_TwrEdge; 
    
    //First Track Specifics
    double Tkr_1_Chisq;
    double Tkr_1_1stChisq;
    double Tkr_1_Gaps;
    double Tkr_1_1stGaps; 
    double Tkr_1_Hits;
    double Tkr_1_1stHits;
    double Tkr_1_1stLayer; 
    
    double Tkr_1_Qual;
    double Tkr_1_Type;
    double Tkr_1_DifHits;
    double Tkr_1_KalEne;
    double Tkr_1_ConEne;
    double Tkr_1_KalThetaMS;
    double Tkr_1_TwrEdge;
    double Tkr_1_PrjTwrEdge;
    double Tkr_1_DieEdge;
    double Tkr_1_xdir;
    double Tkr_1_ydir;
    double Tkr_1_zdir;
    double Tkr_1_Phi;
    double Tkr_1_x0;
    double Tkr_1_y0;
    double Tkr_1_z0;
    
    //Second Track Specifics
    double Tkr_2_Chisq;
    double Tkr_2_1stChisq;
    double Tkr_2_1stGaps; 
    double Tkr_2_Qual;
    double Tkr_2_Type;
    double Tkr_2_Hits;
    double Tkr_2_1stHits;
    double Tkr_2_1stLayer; 
    
    
    double Tkr_2_Gaps;
    double Tkr_2_DifHits;
    double Tkr_2_KalEne;
    double Tkr_2_ConEne;
    double Tkr_2_KalThetaMS;
    double Tkr_2_TwrEdge;
    double Tkr_2_PrjTwrEdge;
    double Tkr_2_DieEdge;
    double Tkr_2_xdir;
    double Tkr_2_ydir;
    double Tkr_2_zdir;
    double Tkr_2_Phi;
    double Tkr_2_x0;
    double Tkr_2_y0;
    double Tkr_2_z0;
};

// Static factory for instantiation of algtool objects
static ToolFactory<TkrValsTool> s_factory;
const IToolFactory& TkrValsToolFactory = s_factory;

// Standard Constructor
TkrValsTool::TkrValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode TkrValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    // get the services
    
    if( serviceLocator() ) {

        // all the XxxValsTools must retrieve EventDataSvc and pass it 
        // to ValBase
        sc = serviceLocator()->service( "EventDataSvc", m_pEventSvc, true );
        if(sc.isFailure()){
            log << MSG::ERROR << "Could not find EventSvc" << endreq;
            return sc;
        }
        setEventSvc(m_pEventSvc);

        sc = serviceLocator()->service( "TkrGeometrySvc", pTkrGeoSvc, true );
        if(sc.isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return sc;
        }
        
        log << MSG::DEBUG << "EventSvc pointer " << m_pEventSvc << endreq;
 
        sc = toolSvc()->retrieveTool("TkrQueryClustersTool", pQueryClusters);
        if (sc.isFailure()) {
            log << MSG::ERROR << "Couldn't retrieve TkrQueryClusterTool" << endreq;
            return StatusCode::FAILURE;
        }
                
        // Which propagator to use?
        int m_PropagatorType = 0; 
        IPropagatorTool* propTool = 0;
        if (m_PropagatorType == 0)
        {
            // Look for the G4PropagatorSvc service
            sc = toolSvc()->retrieveTool("G4PropagatorTool", propTool);
            log << MSG::INFO << "Using Geant4 Particle Propagator" << endreq;
        }
        else
        {
            // Look for GismoGenerator Service
            sc = toolSvc()->retrieveTool("RecoTool", propTool);
            log << MSG::INFO << "Using Gismo Particle Propagator" << endreq;
        }
        pKalParticle = propTool->getPropagator();      
    }
 
    zeroVals();

    // load up the map

    m_ntupleMap["TKR_No_Tracks"]    = &Tkr_No_Tracks;
    m_ntupleMap["TKR_Sum_KalEne"]   = &Tkr_Sum_KalEne;
    m_ntupleMap["TKR_Sum_ConEne"]   = &Tkr_Sum_ConEne;
    m_ntupleMap["TKR_Energy"]       = &Tkr_Energy;
    m_ntupleMap["TKR_Energy_Sum"]   = &Tkr_Energy_Sum;
    m_ntupleMap["TKR_Energy_Corr"]  = &Tkr_Energy_Corr;
    m_ntupleMap["TKR_Edge_Corr"]    = &Tkr_Edge_Corr;
    m_ntupleMap["TKR_Total_Hits"]   = &Tkr_Total_Hits;
    m_ntupleMap["TKR_Thin_Hits"]    = &Tkr_Thin_Hits;
    m_ntupleMap["TKR_Thick_Hits"]   = &Tkr_Thick_Hits;
    m_ntupleMap["TKR_Blank_Hits"]   = &Tkr_Blank_Hits;
    
    m_ntupleMap["TKR_RadLength"]    = &Tkr_RadLength;
    m_ntupleMap["TKR_TwrEdge"]      = &Tkr_TwrEdge;
    
    m_ntupleMap["TKR_1_Chisq"]      = &Tkr_1_Chisq;
    m_ntupleMap["TKR_1_1stChisq"]   = &Tkr_1_1stChisq;
    
    m_ntupleMap["TKR_1_Hits"]       = &Tkr_1_Hits;
    m_ntupleMap["TKR_1_1stHits"]    = &Tkr_1_1stHits;
    m_ntupleMap["TKR_1_1stLayer"]   = &Tkr_1_1stLayer;
    m_ntupleMap["TKR_1_DifHits"]    = &Tkr_1_DifHits;
    
    m_ntupleMap["TKR_1_Gaps"]       = &Tkr_1_Gaps;
    m_ntupleMap["TKR_1_1stGaps"]    = &Tkr_1_1stGaps;
    
    m_ntupleMap["TKR_1_Qual"]       = &Tkr_1_Qual;
    m_ntupleMap["TKR_1_Type"]       = &Tkr_1_Type;
    m_ntupleMap["TKR_1_TwrEdge"]    = &Tkr_1_TwrEdge;
    m_ntupleMap["TKR_1_PrjTwrEdge"] = &Tkr_1_PrjTwrEdge;
    m_ntupleMap["TKR_1_DieEdge"]    = &Tkr_1_DieEdge;
    
    m_ntupleMap["TKR_1_KalEne"]     = &Tkr_1_KalEne;
    m_ntupleMap["TKR_1_ConEne"]     = &Tkr_1_ConEne;
    m_ntupleMap["TKR_1_KalThetaMS"] = &Tkr_1_KalThetaMS;
    
    m_ntupleMap["TKR_1_xdir"]       = &Tkr_1_xdir;
    m_ntupleMap["TKR_1_ydir"]       = &Tkr_1_ydir;
    m_ntupleMap["TKR_1_zdir"]       = &Tkr_1_zdir;
    m_ntupleMap["TKR_1_Phi"]        = &Tkr_1_Phi;
    m_ntupleMap["TKR_1_x0"]         = &Tkr_1_x0;
    m_ntupleMap["TKR_1_y0"]         = &Tkr_1_y0;
    m_ntupleMap["TKR_1_z0"]         = &Tkr_1_z0;
    
    m_ntupleMap["TKR_2_Chisq"]      = &Tkr_2_Chisq;
    m_ntupleMap["TKR_2_1stChisq"]   = &Tkr_2_1stChisq;
    
    m_ntupleMap["TKR_2_Hits"]       = &Tkr_2_Hits;
    m_ntupleMap["TKR_2_1stHits"]    = &Tkr_2_1stHits;
    m_ntupleMap["TKR_2_1stLayer"]   = &Tkr_2_1stLayer;
    m_ntupleMap["TKR_2_DifHits"]    = &Tkr_2_DifHits;
    
    m_ntupleMap["TKR_2_Gaps"]       = &Tkr_2_Gaps;
    m_ntupleMap["TKR_2_1stGaps"]    = &Tkr_2_1stGaps;
    
    m_ntupleMap["TKR_2_Qual"]       = &Tkr_2_Qual;
    m_ntupleMap["TKR_2_Type"]       = &Tkr_2_Type;
    m_ntupleMap["TKR_2_TwrEdge"]    = &Tkr_2_TwrEdge;
    m_ntupleMap["TKR_2_PrjTwrEdge"] = &Tkr_2_PrjTwrEdge;
    m_ntupleMap["TKR_2_DieEdge"]    = &Tkr_2_DieEdge;
    
    m_ntupleMap["TKR_2_KalEne"]     = &Tkr_2_KalEne;
    m_ntupleMap["TKR_2_ConEne"]     = &Tkr_2_ConEne;
    m_ntupleMap["TKR_2_KalThetaMS"] = &Tkr_2_KalThetaMS;
    
    m_ntupleMap["TKR_2_xdir"]       = &Tkr_2_xdir;
    m_ntupleMap["TKR_2_ydir"]       = &Tkr_2_ydir;
    m_ntupleMap["TKR_2_zdir"]       = &Tkr_2_zdir;
    m_ntupleMap["TKR_2_Phi"]        = &Tkr_2_Phi;
    m_ntupleMap["TKR_2_x0"]         = &Tkr_2_x0;
    m_ntupleMap["TKR_2_y0"]         = &Tkr_2_y0;
    m_ntupleMap["TKR_2_z0"]         = &Tkr_2_z0;    
    
    
    return sc;
}

namespace {

    double radThin   = .03; 
    double radThick  = .18; 
    double radTray   = .015;

    // coefs from Miner
    double cfThin    = 0.722;
    double cfThick   = 1.864;
    double cfNoConv  = 0.117;
    double cfRadLen  = 13.07;
    double cfZ       = -0.021;

    double cfNoConv1 = 2.5;

    double max_corr  = 3.0; 
    double rm_hard   = 30.; 
    double rm_soft   = 130;
    double gap       = 18.; 
    double hard_frac = .6; 

    double minHeight = 26.5;
}

StatusCode TkrValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
            
    //Recover EventHeader Pointer
    //SmartDataPtr<Event::EventHeader> pEvent(m_pEventSvc, EventModel::EventHeader);
        
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    SmartDataPtr<Event::TkrVertexCol>     pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    SmartDataPtr<Event::TkrClusterCol> pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    // all variable values are preset to zero. Be sure to re-initialize the ones you care about  

    //Make sure we have valid reconstructed data

    // used to be 89.61; should have been 89.7 // This is the SSD die pitch (width + gap) 
    double die_width = pTkrGeoSvc->ladderNStrips()*pTkrGeoSvc->siStripPitch() +
        + 2.*pTkrGeoSvc->siDeadDistance() + pTkrGeoSvc->ladderGap();

    int nNoConv = pTkrGeoSvc->numNoConverter();
    int nThick  = pTkrGeoSvc->numSuperGlast();
    int nThin   = pTkrGeoSvc->numLayers() - nThick - nNoConv;

    double towerPitch = pTkrGeoSvc->towerPitch(); 

    if (pTracks)
    {   
        // Count number of tracks
        int nParticles = pTracks->size();
        Tkr_No_Tracks   = nParticles;
        
        if(nParticles < 1) return sc;
        
        // Get the first Track - it should be the "Best Track"
        Event::TkrFitConPtr pTrack1 = pTracks->begin();
        
        const Event::TkrFitTrackBase* trackBase = *pTrack1;
        const Event::TkrKalFitTrack* track_1 = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase);
        
        Tkr_1_Chisq        = track_1->getChiSquare();
        Tkr_1_1stChisq     = track_1->chiSquareSegment();
        Tkr_1_1stGaps      = track_1->getNumXFirstGaps() + track_1->getNumYFirstGaps();
        Tkr_1_Qual         = track_1->getQuality();
        Tkr_1_Type         = track_1->getType();
        Tkr_1_Hits         = track_1->getNumHits();
        Tkr_1_1stHits      = track_1->getNumSegmentPoints();
        Tkr_1_1stLayer     = track_1->getLayer();
        Tkr_1_Gaps         = track_1->getNumGaps();
        Tkr_1_KalEne       = track_1->getKalEnergy(); 
        Tkr_1_ConEne       = track_1->getEnergy(); 
        Tkr_1_KalThetaMS   = track_1->getKalThetaMS(); 
        Tkr_1_DifHits      = track_1->getNumXHits()-track_1->getNumYHits();
        
        Point  x1 = track_1->getPosition();
        Vector t1 = track_1->getDirection();
        
        Tkr_1_xdir        = t1.x();
        Tkr_1_ydir        = t1.y();
        Tkr_1_zdir        = t1.z();
        Tkr_1_Phi         = atan(-t1.y()/t1.x()); 
        Tkr_1_x0          = x1.x();
        Tkr_1_y0          = x1.y();
        Tkr_1_z0          = x1.z();
        
        double z_dist    = fabs((pTkrGeoSvc->trayHeight()+3.)/t1.z()); 
        double x_twr = sign(x1.x())*(fmod(fabs(x1.x()),towerPitch) - towerPitch/2.);
        double y_twr = sign(x1.y())*(fmod(fabs(x1.y()),towerPitch) - towerPitch/2.);
        double x_prj = x_twr - t1.x()*z_dist;
        double y_prj = y_twr - t1.y()*z_dist; 
        
        Tkr_1_TwrEdge    = (fabs(x_twr) > fabs(y_twr)) ? fabs(x_twr) : fabs(y_twr);
        Tkr_1_PrjTwrEdge = (fabs(x_prj) > fabs(y_prj)) ? fabs(x_prj) : fabs(y_prj);
        Tkr_1_TwrEdge    = towerPitch/2. - Tkr_1_TwrEdge;
        Tkr_1_PrjTwrEdge = towerPitch/2. - Tkr_1_PrjTwrEdge;
        
        double x_die = sign(x_twr)*(fmod(fabs(x_twr),die_width) - die_width/2.);
        double y_die = sign(y_twr)*(fmod(fabs(y_twr),die_width) - die_width/2.);
        
        Tkr_1_DieEdge  = (fabs(x_die) > fabs(y_die)) ? fabs(x_die) : fabs(y_die);
        Tkr_1_DieEdge  = die_width/2. - Tkr_1_DieEdge; 
        
        if(nParticles > 1) {
            pTrack1++;
            trackBase = *pTrack1;
            const Event::TkrKalFitTrack* track_2 = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase);
            
            Tkr_2_Chisq        = track_2->getChiSquare();
            Tkr_2_1stChisq     = track_2->chiSquareSegment();
            Tkr_2_1stGaps      = track_2->getNumXFirstGaps() + track_2->getNumYFirstGaps();
            Tkr_2_Qual         = track_2->getQuality();
            Tkr_2_Type         = track_2->getType();
            Tkr_2_Hits         = track_2->getNumHits();
            Tkr_2_1stHits      = track_2->getNumSegmentPoints();
            Tkr_2_1stLayer     = track_2->getLayer();
            Tkr_2_Gaps         = track_2->getNumGaps();
            Tkr_2_KalEne       = track_2->getKalEnergy(); 
            Tkr_2_ConEne       = track_2->getEnergy(); 
            Tkr_2_KalThetaMS   = track_2->getKalThetaMS(); 
            Tkr_2_DifHits      = track_2->getNumXHits()-track_2->getNumYHits();
            
            Point  x2 = track_2->getPosition();
            Vector t2 = track_2->getDirection();
            Tkr_2_xdir       = t2.x();
            Tkr_2_ydir       = t2.y();
            Tkr_2_zdir       = t2.z();
            Tkr_2_Phi        = atan(-t2.y()/t2.x()); 
            Tkr_2_x0         = x2.x();
            Tkr_2_y0         = x2.y();
            Tkr_2_z0         = x2.z();
            
            double x_twr = sign(x2.x())*(fmod(fabs(x2.x()),towerPitch) - towerPitch/2.);
            double y_twr = sign(x2.y())*(fmod(fabs(x2.y()),towerPitch) - towerPitch/2.);
            double x_prj = x_twr - t2.x()*z_dist;
            double y_prj = y_twr - t2.y()*z_dist; 
            
            Tkr_2_TwrEdge    = (fabs(x_twr) > fabs(y_twr)) ? fabs(x_twr) : fabs(y_twr);
            Tkr_2_PrjTwrEdge = (fabs(x_prj) > fabs(y_prj)) ? fabs(x_prj) : fabs(y_prj);
            Tkr_2_TwrEdge    = towerPitch/2. - Tkr_2_TwrEdge;
            Tkr_2_PrjTwrEdge = towerPitch/2. - Tkr_2_PrjTwrEdge;
            
            double x_die = sign(x_twr)*(fmod(fabs(x_twr),die_width) - die_width/2.);
            double y_die = sign(y_twr)*(fmod(fabs(y_twr),die_width) - die_width/2.);
            
            Tkr_2_DieEdge  = (fabs(x_die) > fabs(y_die)) ? fabs(x_die) : fabs(y_die);
            Tkr_2_DieEdge  = die_width/2. - Tkr_2_DieEdge; 
        }
        
        Tkr_Sum_KalEne    = Tkr_1_KalEne+Tkr_2_KalEne; 
        Tkr_Sum_ConEne    = Tkr_1_ConEne+Tkr_2_ConEne;      
        
        
        // Computation of the tracker contribution to the total energy 
        double costh = fabs(t1.z()); 
        double arc_min = (x1.z() + minHeight)/costh; 
        pKalParticle->setStepStart(x1, t1, arc_min);
        double total_radlen = pKalParticle->radLength(); 
        
        // Compute the sum-of radiation_lengths x Hits in each layer
        double tracker_ene_sum  = 0.;
        double tracker_ene_corr = 0.; 
        double rad_len_sum  = 0.; 
        double radlen       = 0.;
        double radlen_old   = 0.; 
        double arc_len      = 0.; 
        int    total_hits   = 0; 
        int    thin_hits    = 0;
        int    thick_hits   = 0; 
        int    blank_hits   = 0; 
        double ave_edge     = 0.; 
        
        int top_plane     = track_1->getLayer(); 
        
        int max_planes = pTkrGeoSvc->numLayers();
        
        double two_phi = 2.*Tkr_1_Phi; 
        double angle_factor_G = 1.- (1.- costh)*sin(two_phi)*sin(two_phi); 
        double angle_factor_R = 1.+ (1./costh - 1.)*sqrt(0.5)*sin(two_phi)*sin(two_phi); 
        gap     *= angle_factor_G;
        rm_hard *= angle_factor_R;
        rm_soft *= angle_factor_R;

        // take care of conversion 1/2 way through first radiator
        radlen = radThin/2.;
        if(top_plane >= nThin) {
            // "0." won't happen for standard geometry, because 3 layers are required for track
            radlen = (top_plane<max_planes-nNoConv ? 0.5*radThick : 0.);
        }
        
        for(int iplane = top_plane; iplane < max_planes; iplane++) {
            
            double xms = 0.;
            double yms = 0.;
                        
            if(iplane > top_plane) {
                Event::TkrFitMatrix Q = pKalParticle->mScat_Covr(Tkr_Sum_ConEne/2., arc_len);
                xms = Q.getcovX0X0();
                yms = Q.getcovY0Y0();
                radlen += pKalParticle->radLength(arc_len); 
            }
            double xSprd = sqrt(4.+xms*16.); // 4.0 sigma and not smaller then 2mm (was 2.5 sigma)
            double ySprd = sqrt(4.+yms*16.); // Limit to a tower... 
            if(xSprd > pTkrGeoSvc->trayWidth()/2.) xSprd = pTkrGeoSvc->trayWidth()/2.;
            if(ySprd > pTkrGeoSvc->trayWidth()/2.) ySprd = pTkrGeoSvc->trayWidth()/2.;
            
            // Assume location of shower center is given by 1st track
            Point x_hit = x1 + arc_len*t1;
            int numHits = pQueryClusters->numberOfHitsNear(iplane, xSprd, ySprd, x_hit);
            
            int outside = 0; 
            int iView   = 0;
            double layer_edge = twrEdgeT(x_hit.x(), x_hit.y(), towerPitch, iView, outside);
            double rm_frac_plus = (layer_edge-gap/2.)/rm_soft; 
            double in_frac_soft = circle_fracT(rm_frac_plus);
            if(!outside) {
                double rm_frac_minus = (layer_edge + gap/2.)/rm_soft;
                if(rm_frac_minus > 0.)in_frac_soft += circle_fracT(-rm_frac_minus);
            }
            if(in_frac_soft < .01) in_frac_soft = .01; 
            rm_frac_plus = (layer_edge-gap/2.)/rm_hard; 
            double in_frac_hard = circle_fracT(rm_frac_plus);
            if(!outside) {
                double rm_frac_minus = (layer_edge + gap/2.)/rm_hard;
                if(rm_frac_minus > 0.)in_frac_hard += circle_fracT(-rm_frac_minus);
            }
            if(in_frac_hard < .01) in_frac_hard = .01; 
            double corr_factor = 1./((1.-hard_frac)*in_frac_soft + hard_frac*in_frac_hard);
            if(corr_factor > max_corr) corr_factor = max_corr; 
            double delta_rad= radlen-radlen_old;
            if(iplane < nThin) {
                if(delta_rad < radThin/costh) 
                    delta_rad=(radThin+radTray)/costh;
            }
            else if(iplane < max_planes-2) {
                if(delta_rad < radThick/costh) 
                    delta_rad=(radThick+radTray)/costh;
            }
            if(iplane < nThin)                   thin_hits  += numHits;
            else if(iplane < max_planes-nNoConv) thick_hits += numHits;
            else                                 blank_hits += numHits;
            
            double ene      = delta_rad*numHits*10.;
            double ene_corr = corr_factor*delta_rad; //*ene;
            tracker_ene_sum  += ene; 
            tracker_ene_corr += ene_corr;
            total_hits       += numHits; 
            ave_edge         += layer_edge*delta_rad; 
            rad_len_sum      += delta_rad;
            
            // Increment arc-length
            //arc_len += pTkrGeoSvc->trayHeight()/fabs(dir_ini.z());
            int nextPlane = iplane+1;
            if (iplane==max_planes-1) nextPlane--;
            double deltaZ = pTkrGeoSvc->getReconLayerZ(iplane) -
                pTkrGeoSvc->getReconLayerZ(nextPlane);
            arc_len += fabs( deltaZ/t1.z()); 
            radlen_old = radlen; 
        }
        Tkr_RadLength  = rad_len_sum;
        // Coef's from Lin. Regression analysis in Miner
        // defined in anonymous namespace above
        Tkr_Energy     = (cfThin*thin_hits + cfThick*thick_hits + cfNoConv*blank_hits
            + cfRadLen*rad_len_sum + cfZ*x1.z())/costh;
        Tkr_Energy_Sum = tracker_ene_sum + cfNoConv1*blank_hits;  
        Tkr_Edge_Corr  = tracker_ene_corr/rad_len_sum;
        Tkr_Energy_Corr= Tkr_Edge_Corr*Tkr_Energy;
        Tkr_Total_Hits = total_hits;
        Tkr_Thin_Hits  = thin_hits;
        Tkr_Thick_Hits = thick_hits;
        Tkr_Blank_Hits = blank_hits; 
        Tkr_TwrEdge    = ave_edge/rad_len_sum; 
    }          
    
    return sc;
}
