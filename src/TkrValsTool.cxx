/** @file TkrValsTool.cxx
@brief Calculates the Tkr analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrValsTool.cxx,v 1.19 2003/06/12 00:00:06 atwood Exp $
*/

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
#include "Event/Recon/TkrRecon/TkrFitMatrix.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
//#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GaudiKernel/IToolSvc.h"

#ifndef M_PI
#define M_PI = 3.14159265358979323846
#endif

namespace {
    const int _nTowers = 4;
    
    double sign(double x) { return (x<0.) ? -1.:1.;} 
    
    double twrEdgeT(double x, double y, int nXTowers, int nYTowers, 
        double pitch, int &XY, int &outer) {
        double edge = 0.; 
        double x_twr = sign(x)*(fmod(fabs(x),pitch) - pitch/2.);
        double y_twr = sign(y)*(fmod(fabs(y),pitch) - pitch/2.);
        
        outer = 0; 
        int outerX = 0;
        int outerY = 0;
        
        if(fabs(x_twr) > fabs(y_twr)) {
            edge = pitch/2. - fabs(x_twr);
            XY = 1; 
            if(fabs(x) > 0.5*(nXTowers-1)*pitch) outer = 1;
        }
        else {
            edge = pitch/2. - fabs(y_twr);
            XY = 2;
            if(fabs(y) > 0.5*(nYTowers-1)*pitch) outer = 1;
        }
        return edge;
    }
    
	double circle_frac(double r) {
        double rl = (fabs(r) < 1.) ? fabs(r):1.; 
        double a_slice = 2.*(M_PI/4. - rl*sqrt(1.-rl*rl)/2. - asin(rl)/2.);
        double in_frac = 1.-a_slice/M_PI;
        if(r < 0.) in_frac = a_slice/M_PI;
        return in_frac;
    }
    
    double circle_frac_simp(double r, double angle_factor) {
        double slice_0 = circle_frac(r);
        double slice_p = circle_frac(r+angle_factor);
        double slice_m = circle_frac(r-angle_factor);
        return (slice_p + 4.*slice_0 + slice_m)/6.;
    }
    
    double contained_frac(double x, double y, double pitch, double gap,  
        double r, double costh, double phi) {
        // Get the projected angles for the gap
        double tanth = sqrt(1.-costh*costh)/costh;
        double gap_x = gap - 35.*sin(phi)*tanth; 
        if(gap_x < 5.) gap_x = 5.;
        double gap_y = gap - 35.*cos(phi)*tanth;
        if(gap_y < 5.) gap_y = 5.; 
        
        // X Edges
        double x_twr = sign(x)*(fmod(fabs(x),pitch) - pitch/2.);
        double edge = pitch/2. - fabs(x_twr);
        double r_frac_plus = (edge-gap_x/2.)/r; 
        double angle_factor = sin(phi)*(1./costh - 1.);
        double in_frac_x  =  circle_frac_simp(r_frac_plus, angle_factor);
        if(fabs(x) < 1.5*pitch) { // X edge is not outside limit of LAT
            double r_frac_minus = (edge + gap_x/2.)/r;
            in_frac_x += circle_frac_simp(-r_frac_minus, angle_factor);
        }
        
        // Y Edges
        double y_twr = sign(y)*(fmod(fabs(y),pitch) - pitch/2.);
        edge = pitch/2. - fabs(y_twr);
        r_frac_plus = (edge-gap_y/2.)/r; 
        angle_factor = cos(phi)*(1./costh - 1.);
        double in_frac_y  =  circle_frac_simp(r_frac_plus, angle_factor);
        if(fabs(y) < 1.5*pitch) { // X edge is not outside limit of LAT
            double r_frac_minus = (edge + gap_y/2.)/r;
            in_frac_y += circle_frac_simp(-r_frac_minus, angle_factor);
        }
        
        // Cross term assumes x and y are independent 
        double in_frac = 1.;
        if(in_frac_x > .999) in_frac = in_frac_y;
        else if(in_frac_y > .999) in_frac = in_frac_x; 
        else in_frac = 1. - (1.-in_frac_x) - (1.-in_frac_y) + 
            (1.-in_frac_x)*(1.-in_frac_y);  //Cross Term Correction?  
        if(in_frac < .01) in_frac = .01;
        
        return in_frac;
    }
}

/*! @class TkrValsTool
@brief calculates Tkr values

@authors Bill Atwood, Leon Rochester
*/

class TkrValsTool :  public ValBase 
{
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
    /// pointer to queryclusterstool
    ITkrQueryClustersTool* pQueryClusters;
    /// 
    IPropagatorSvc* m_propSvc;
    
    IKalmanParticle*       pKalParticle;
    
    
    //Global Track Tuple Items
    double Tkr_No_Tracks;
    double Tkr_Sum_KalEne; 
    double Tkr_Sum_ConEne;
    double Tkr_Energy;
    double Tkr_Energy_Sum;
    double Tkr_Energy_Corr;
    double Tkr_Edge_Corr;
    double Tkr_HDCount; 
    double Tkr_Total_Hits;
    double Tkr_Thin_Hits;
    double Tkr_Thick_Hits;
    double Tkr_Blank_Hits;  
    double Tkr_RadLength; 
    double Tkr_TwrEdge; 
    double Tkr_TrackLength;
    
    //First Track Specifics
    double Tkr_1_Chisq;
    double Tkr_1_FirstChisq;
    double Tkr_1_Gaps;
    double Tkr_1_FirstGaps; 
    double Tkr_1_Hits;
    double Tkr_1_FirstHits;
    double Tkr_1_FirstLayer; 
    
    double Tkr_1_Qual;
    double Tkr_1_Type;

    double Tkr_1_DifHits;
    double Tkr_1_KalEne;
    double Tkr_1_ConEne;
    double Tkr_1_KalThetaMS;
    double Tkr_1_TwrEdge;
    double Tkr_1_PrjTwrEdge;
    double Tkr_1_DieEdge;
	double Tkr_1_TwrGap;
    double Tkr_1_xdir;
    double Tkr_1_ydir;
    double Tkr_1_zdir;
    double Tkr_1_Phi;
	double Tkr_1_Theta;
    double Tkr_1_x0;
    double Tkr_1_y0;
    double Tkr_1_z0;
	double Tkr_1_Sxx;
	double Tkr_1_Sxy;
	double Tkr_1_Syy;
	double Tkr_1_DTheta;
	double Tkr_1_DPhi;
    
    //Second Track Specifics
    double Tkr_2_Chisq;
    double Tkr_2_FirstChisq;
    double Tkr_2_FirstGaps; 
    double Tkr_2_Qual;
    double Tkr_2_Type;
    double Tkr_2_Hits;
    double Tkr_2_FirstHits;
    double Tkr_2_FirstLayer; 
    
    
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
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode TkrValsTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;
    
    MsgStream log(msgSvc(), name());
    
    if((ValBase::initialize()).isFailure()) return StatusCode::FAILURE;
    
    // get the services
    
    if( serviceLocator() ) {
        
        if(service( "TkrGeometrySvc", pTkrGeoSvc, true ).isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return fail;
        }
               
        // pick up the chosen propagator
        if (service("GlastPropagatorSvc", m_propSvc, true ).isFailure()) {
            log << MSG::ERROR << "Couldn't find the GlastPropagatorSvc!" << endreq;
            return fail;
        }
        
        pKalParticle = m_propSvc->getPropagator();
        
    } else {
        return fail;
    }
    
    if (toolSvc()->retrieveTool("TkrQueryClustersTool", pQueryClusters).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrQueryClusterTool" << endreq;
        return fail;
    }

    // load up the map
    
    addItem("TkrNumTracks",   &Tkr_No_Tracks);
    addItem("TkrSumKalEne",   &Tkr_Sum_KalEne);
    addItem("TkrSumConEne",   &Tkr_Sum_ConEne);
    addItem("TkrEnergy",      &Tkr_Energy);
    addItem("TkrEnergySum",   &Tkr_Energy_Sum);
    addItem("TkrEnergyCorr",  &Tkr_Energy_Corr);
    addItem("TkrEdgeCorr",    &Tkr_Edge_Corr);
    addItem("TkrHDCount",     &Tkr_HDCount); 
    addItem("TkrTotalHits",   &Tkr_Total_Hits);
    addItem("TkrThinHits",    &Tkr_Thin_Hits);
    addItem("TkrThickHits",   &Tkr_Thick_Hits);
    addItem("TkrBlankHits",   &Tkr_Blank_Hits);
    
    addItem("TkrRadLength",   &Tkr_RadLength);
    addItem("TkrTwrEdge",     &Tkr_TwrEdge);
    addItem("TkrTrackLength", &Tkr_TrackLength);
    
    addItem("Tkr1Chisq",      &Tkr_1_Chisq);
    addItem("Tkr1FirstChisq", &Tkr_1_FirstChisq);
    addItem("Tkr1Hits",       &Tkr_1_Hits);
    addItem("Tkr1FirstHits",  &Tkr_1_FirstHits);
    addItem("Tkr1FirstLayer", &Tkr_1_FirstLayer);
    addItem("Tkr1DifHits",    &Tkr_1_DifHits);
    
    addItem("Tkr1Gaps",       &Tkr_1_Gaps);
    addItem("Tkr1FirstGaps",  &Tkr_1_FirstGaps);
    
    addItem("Tkr1Qual",       &Tkr_1_Qual);
    addItem("Tkr1Type",       &Tkr_1_Type);
    addItem("Tkr1TwrEdge",    &Tkr_1_TwrEdge);
    addItem("Tkr1PrjTwrEdge", &Tkr_1_PrjTwrEdge);
    addItem("Tkr1DieEdge",    &Tkr_1_DieEdge);
    addItem("Tkr1TwrGap",     &Tkr_1_TwrGap);
    
    addItem("Tkr1KalEne",     &Tkr_1_KalEne);
    addItem("Tkr1ConEne",     &Tkr_1_ConEne);
    addItem("Tkr1KalThetaMS", &Tkr_1_KalThetaMS);
    
    addItem("Tkr1XDir",       &Tkr_1_xdir);
    addItem("Tkr1YDir",       &Tkr_1_ydir);
    addItem("Tkr1ZDir",       &Tkr_1_zdir);
    addItem("Tkr1Phi",        &Tkr_1_Phi);
	addItem("Tkr1Theta",      &Tkr_1_Theta);
    addItem("Tkr1X0",         &Tkr_1_x0);
    addItem("Tkr1Y0",         &Tkr_1_y0);
    addItem("Tkr1Z0",         &Tkr_1_z0);

	addItem("Tkr1DTheta",     &Tkr_1_DTheta);
    addItem("Tkr1DPhi",       &Tkr_1_DPhi);
	addItem("Tkr1SXX",        &Tkr_1_Sxx);
    addItem("Tkr1SXY",        &Tkr_1_Sxy);
    addItem("Tkr1SYY",        &Tkr_1_Syy);

    addItem("Tkr2Chisq",      &Tkr_2_Chisq);
    addItem("Tkr2FirstChisq", &Tkr_2_FirstChisq);
    
    addItem("Tkr2Hits",       &Tkr_2_Hits);
    addItem("Tkr2FirstHits",  &Tkr_2_FirstHits);
    addItem("Tkr2FirstLayer", &Tkr_2_FirstLayer);
    addItem("Tkr2DifHits",    &Tkr_2_DifHits);
    
    addItem("Tkr2Gaps",       &Tkr_2_Gaps);
    addItem("Tkr2FirstGaps",  &Tkr_2_FirstGaps);
    
    addItem("Tkr2Qual",       &Tkr_2_Qual);
    addItem("Tkr2Type",       &Tkr_2_Type);
    addItem("Tkr2TwrEdge",    &Tkr_2_TwrEdge);
    addItem("Tkr2PrjTwrEdge", &Tkr_2_PrjTwrEdge);
    addItem("Tkr2DieEdge",    &Tkr_2_DieEdge);
   
    addItem("Tkr2KalEne",     &Tkr_2_KalEne);
    addItem("Tkr2ConEne",     &Tkr_2_ConEne);
    addItem("Tkr2KalThetaMS", &Tkr_2_KalThetaMS);
    
    addItem("Tkr2XDir",       &Tkr_2_xdir);
    addItem("Tkr2YDir",       &Tkr_2_ydir);
    addItem("Tkr2ZDir",       &Tkr_2_zdir);
    addItem("Tkr2Phi",        &Tkr_2_Phi);
    addItem("Tkr2X0",         &Tkr_2_x0);
    addItem("Tkr2Y0",         &Tkr_2_y0);
    addItem("Tkr2Z0",         &Tkr_2_z0);    
    
    zeroVals();
    
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
    double hard_frac = .7; 
}

StatusCode TkrValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    //Recover EventHeader Pointer
    //SmartDataPtr<Event::EventHeader> pEvent(m_pEventSvc, EventModel::EventHeader);
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    //SmartDataPtr<Event::TkrVertexCol>     pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    //SmartDataPtr<Event::TkrClusterCol> pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    if(!pTracks) return StatusCode::FAILURE;

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
        Tkr_1_FirstChisq   = track_1->chiSquareSegment();
        Tkr_1_FirstGaps    = track_1->getNumXFirstGaps() + track_1->getNumYFirstGaps();
        Tkr_1_Qual         = track_1->getQuality();
        Tkr_1_Type         = track_1->getType();
        Tkr_1_Hits         = track_1->getNumHits();
        Tkr_1_FirstHits    = track_1->getNumSegmentPoints();
        Tkr_1_FirstLayer   = track_1->getLayer();
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

        Tkr_1_x0          = x1.x();
        Tkr_1_y0          = x1.y();
        Tkr_1_z0          = x1.z();

		Tkr_1_Phi         = atan(-t1.y()/t1.x()); 
		Tkr_1_Theta       = acos(t1.z());

		Event::TkrFitMatrix  Tkr_1_Cov = track_1->getTrackCov();
		Tkr_1_Sxx         = Tkr_1_Cov.getcovSxSx();
		Tkr_1_Sxy         = Tkr_1_Cov.getcovSxSy();
		Tkr_1_Syy         = Tkr_1_Cov.getcovSySy();
		double tzmax      = (fabs(t1.z()) > .999999)? fabs(t1.z()) : .999999;
		double tzTerm     = tzmax/(1.-tzmax*tzmax);
		Tkr_1_DTheta      = sqrt((t1.x()*t1.z()*t1.z())*(t1.x()*t1.z()*t1.z())*Tkr_1_Sxx +
                                 (t1.y()*t1.z()*t1.z())*(t1.y()*t1.z()*t1.z())*Tkr_1_Syy);
                               //  + (t1.x()*t1.z()*t1.z())*(t1.y()*t1.z()*t1.z())*Tkr1_1_Sxy);
		Tkr_1_DPhi        = sqrt((t1.x()*tzTerm)*(t1.x()*tzTerm)*Tkr_1_Sxx +
                                 (t1.y()*tzTerm)*(t1.y()*tzTerm)*Tkr_1_Syy);
		                      //  + (t1.x()*sxsyTerm)*(t1.y()*sxsyTerm)*Tkr1_1_Sxx); 

        Tkr_TrackLength = -Tkr_1_z0/Tkr_1_zdir;
        
        double z_dist    = fabs((pTkrGeoSvc->trayHeight()+3.)/t1.z()); 
        double x_twr = sign(x1.x())*(fmod(fabs(x1.x()),towerPitch) - towerPitch/2.);
        double y_twr = sign(x1.y())*(fmod(fabs(x1.y()),towerPitch) - towerPitch/2.);
        double x_prj = x_twr - t1.x()*z_dist;
        double y_prj = y_twr - t1.y()*z_dist; 
        
        Tkr_1_TwrEdge    = (fabs(x_twr) > fabs(y_twr)) ? fabs(x_twr) : fabs(y_twr);
        Tkr_1_PrjTwrEdge = (fabs(x_prj) > fabs(y_prj)) ? fabs(x_prj) : fabs(y_prj);
        Tkr_1_TwrEdge    = towerPitch/2. - Tkr_1_TwrEdge;
        Tkr_1_PrjTwrEdge = towerPitch/2. - Tkr_1_PrjTwrEdge;
        
		// New section go compute gap lengths in tracker and cal
		double x_slope   = (fabs(t1.x()) > .0001)? t1.x():.00001;
		double s_x       = (sign(t1.x())*towerPitch/2. - x_twr)/x_slope; 
	    double y_slope   = (fabs(t1.y()) > .0001)? t1.y():.00001;
		double s_y       = (sign(t1.y())*towerPitch/2. - y_twr)/y_slope;
		
		Tkr_1_TwrGap = 0.; 
        if(s_x < s_y) { // Goes out x side of CAL Module
		    if(x1.z() + s_x*t1.z() > 0) {
			    double s_max = fabs((x1.z())/t1.z());
			    Tkr_1_TwrGap = gap/fabs(x_slope);
			    if((Tkr_1_TwrGap + s_x)> s_max ) Tkr_1_TwrGap = s_max-s_x;
		    }
	    }
	    else {          // Goes out y side
		    if(x1.z() + s_y*t1.z() > 0) {
			    double s_max = fabs(x1.z()/t1.z());
			    Tkr_1_TwrGap = gap/fabs(y_slope);
			    if((Tkr_1_TwrGap + s_y)> s_max ) Tkr_1_TwrGap = s_max-s_y;
		    }
	    }

		// SSD Die loaction and edge... 
        double x_die = sign(x_twr)*(fmod(fabs(x_twr),die_width) - die_width/2.);
        double y_die = sign(y_twr)*(fmod(fabs(y_twr),die_width) - die_width/2.);
        
        Tkr_1_DieEdge  = (fabs(x_die) > fabs(y_die)) ? fabs(x_die) : fabs(y_die);
        Tkr_1_DieEdge  = die_width/2. - Tkr_1_DieEdge; 
        
        if(nParticles > 1) {
            pTrack1++;
            trackBase = *pTrack1;
            const Event::TkrKalFitTrack* track_2 = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase);
            
            Tkr_2_Chisq        = track_2->getChiSquare();
            Tkr_2_FirstChisq     = track_2->chiSquareSegment();
            Tkr_2_FirstGaps      = track_2->getNumXFirstGaps() + track_2->getNumYFirstGaps();
            Tkr_2_Qual         = track_2->getQuality();
            Tkr_2_Type         = track_2->getType();
            Tkr_2_Hits         = track_2->getNumHits();
            Tkr_2_FirstHits      = track_2->getNumSegmentPoints();
            Tkr_2_FirstLayer     = track_2->getLayer();
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
        double arc_min = (x1.z() - pTkrGeoSvc->calZTop())/costh; 
        pKalParticle->setStepStart(x1, t1, arc_min);
        //double total_radlen = pKalParticle->radLength(); 
        
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
        
        for(int iplane = top_plane; iplane < max_planes; iplane++) {
            
            double xms = 0.;
            double yms = 0.;
            
            if(iplane > top_plane) {
                Event::TkrFitMatrix Q = pKalParticle->mScat_Covr(Tkr_Sum_ConEne/2., arc_len);
                xms = Q.getcovX0X0();
                yms = Q.getcovY0Y0();
                radlen = pKalParticle->radLength(arc_len); 
            }
            double xSprd = sqrt(4.+xms*16.); // 4.0 sigma and not smaller then 2mm (was 2.5 sigma)
            double ySprd = sqrt(4.+yms*16.); // Limit to a tower... 
            double halfTray = 0.5*pTkrGeoSvc->trayWidth();
            if(xSprd > halfTray) xSprd = halfTray;
            if(ySprd > halfTray) ySprd = halfTray   ;
            
            // Assume location of shower center is given by 1st track
            Point x_hit = x1 + arc_len*t1;
            int numHits = pQueryClusters->numberOfHitsNear(iplane, xSprd, ySprd, x_hit);
            if(iplane == top_plane) {
                double xRgn = 30.*sqrt(1+(cos(Tkr_1_Phi)/costh)*(cos(Tkr_1_Phi)/costh));
                double yRgn = 30.*sqrt(1+(sin(Tkr_1_Phi)/costh)*(sin(Tkr_1_Phi)/costh)); 
                Tkr_HDCount = pQueryClusters->numberOfUUHitsNear(iplane, xRgn, yRgn, x_hit);
            }
            
            int outside = 0; 
            int iView   = 0;
            double layer_edge = twrEdgeT(x_hit.x(), x_hit.y(), 
                                         pTkrGeoSvc->numXTowers(), pTkrGeoSvc->numYTowers(), 
                                         towerPitch, iView, outside);
 
            double in_frac_soft = contained_frac(x_hit.x(), x_hit.y(), towerPitch, gap,  
                                                  rm_soft, costh, Tkr_1_Phi);
            if(in_frac_soft < .01) in_frac_soft = .01; 
 
            double in_frac_hard = contained_frac(x_hit.x(), x_hit.y(), towerPitch, gap,  
                                                  rm_hard, costh, Tkr_1_Phi);
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

		// take care of conversion 1/2 way through first radiator
        if(top_plane >= nThin) {
            // "0." won't happen for standard geometry, because 3 layers are required for track
            rad_len_sum += (top_plane<max_planes-nNoConv ? 0.5*radThick : 0.);
        }
		else {
            rad_len_sum += radThin/2.;
		}
        Tkr_RadLength  = rad_len_sum;


    }          
    
    return sc;
}
