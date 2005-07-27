/** @file TkrValsTool.cxx
@brief Calculates the Tkr analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrValsTool.cxx,v 1.59 2005/06/17 05:11:16 atwood Exp $
*/

// To Do:
// implement better code to check if in tower
// Don't forget to remove the "1.5"s!! Done
// xEdge and yEdge... Done

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
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "geometry/Ray.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GaudiKernel/IToolSvc.h"
#include "geometry/Ray.h"


// M_PI defined in ValBase.h

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

    double towerEdge(Point pos) const;
    double containedFraction(Point pos, double gap, double r, double costh, double phi) const;

    // some local constants
    double m_towerPitch;
    int    m_xNum;
    int    m_yNum;

    // some pointers to services

    /// pointer to tracker geometry
    ITkrGeometrySvc*       m_tkrGeom;
    /// GlastDetSvc used for access to detector info
    IGlastDetSvc*         m_detSvc; 
    /// pointer to queryclusterstool
    ITkrQueryClustersTool* pQueryClusters;
    /// 
    IPropagatorSvc* m_propSvc;

    IPropagator*           m_G4PropTool;    

    //Global Track Tuple Items
    float Tkr_No_Tracks;
    float Tkr_Sum_KalEne; 
    float Tkr_Sum_ConEne;
    float Tkr_Energy;
    float Tkr_Energy_Corr;
    float Tkr_HDCount; 
    float Tkr_Total_Hits;
    float Tkr_Thin_Hits;
    float Tkr_Thick_Hits;
    float Tkr_Blank_Hits;  
    float Tkr_RadLength; 
    float Tkr_TwrEdge; 
    float Tkr_TrackLength;

    //First Track Specifics
    float Tkr_1_Chisq;
    float Tkr_1_FirstChisq;
    float Tkr_1_Gaps;
	float Tkr_1_FirstGapPlane; 
	float Tkr_1_GapX;
	float Tkr_1_GapY;
    float Tkr_1_FirstGaps; 
    float Tkr_1_Hits;
    float Tkr_1_FirstHits;
    float Tkr_1_FirstLayer; 
    float Tkr_1_LastLayer; 

    float Tkr_1_Qual;
    float Tkr_1_Type;

    float Tkr_1_DifHits;
    float Tkr_1_KalEne;
    float Tkr_1_ConEne;
    float Tkr_1_KalThetaMS;
    float Tkr_1_TwrEdge;
    float Tkr_1_PrjTwrEdge;
    float Tkr_1_DieEdge;
    float Tkr_1_TwrGap;
    float Tkr_1_xdir;
    float Tkr_1_ydir;
    float Tkr_1_zdir;
    float Tkr_1_Phi;
    float Tkr_1_Theta;
    float Tkr_1_x0;
    float Tkr_1_y0;
    float Tkr_1_z0;
    float Tkr_1_Sxx;
    float Tkr_1_Sxy;
    float Tkr_1_Syy;
    float Tkr_1_ThetaErr;
    float Tkr_1_PhiErr;
    float Tkr_1_ErrAsym;
    float Tkr_1_CovDet;
    float Tkr_1_ToTFirst;
    float Tkr_1_ToTAve;
    float Tkr_1_ToTTrAve;
    float Tkr_1_ToTAsym;
    float Tkr_1_ChisqAsym;
    float Tkr_1_SSDVeto; 

    //Second Track Specifics
    float Tkr_2_Chisq;
    float Tkr_2_FirstChisq;
    float Tkr_2_FirstGaps; 
    float Tkr_2_Qual;
    float Tkr_2_Type;
    float Tkr_2_Hits;
    float Tkr_2_FirstHits;
    float Tkr_2_FirstLayer; 
    float Tkr_2_LastLayer; 

    float Tkr_2_Gaps;
    float Tkr_2_DifHits;
    float Tkr_2_KalEne;
    float Tkr_2_ConEne;
    float Tkr_2_KalThetaMS;
    float Tkr_2_TwrEdge;
    float Tkr_2_PrjTwrEdge;
    float Tkr_2_DieEdge;
    float Tkr_2_xdir;
    float Tkr_2_ydir;
    float Tkr_2_zdir;
    float Tkr_2_Phi;
    float Tkr_2_Theta;
    float Tkr_2_x0;
    float Tkr_2_y0;
    float Tkr_2_z0;

    float Tkr_2TkrAngle;
    float Tkr_2TkrHDoca;

    // here's some test stuff... if it works for a couple it will work for all
    //float Tkr_float;
    //int   Tkr_int;
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

        if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return fail;
        }

        m_towerPitch = m_tkrGeom->towerPitch();
        m_xNum       = m_tkrGeom->numXTowers();
        m_yNum       = m_tkrGeom->numYTowers();

        // find GlastDevSvc service
        if (service("GlastDetSvc", m_detSvc, true).isFailure()){
            log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
            return StatusCode::FAILURE;
        }

        IToolSvc* toolSvc = 0;
        if(service("ToolSvc", toolSvc, true).isFailure()) {
            log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
            return StatusCode::FAILURE;
        }
        if(!toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)) {
            log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
            return StatusCode::FAILURE;
        }

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
    addItem("TkrEnergyCorr",  &Tkr_Energy_Corr);
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
    addItem("Tkr1LastLayer",  &Tkr_1_LastLayer);
    addItem("Tkr1DifHits",    &Tkr_1_DifHits);

    addItem("Tkr1Gaps",       &Tkr_1_Gaps);
    addItem("Tkr1FirstGapPlane",&Tkr_1_FirstGapPlane);
	addItem("Tkr1XGap",       &Tkr_1_GapX);
    addItem("Tkr1YGap",       &Tkr_1_GapY);
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

    addItem("Tkr1ThetaErr",   &Tkr_1_ThetaErr);
    addItem("Tkr1PhiErr",     &Tkr_1_PhiErr);
    addItem("Tkr1ErrAsym",    &Tkr_1_ErrAsym);
    addItem("Tkr1CovDet",     &Tkr_1_CovDet);
    addItem("Tkr1SXX",        &Tkr_1_Sxx);
    addItem("Tkr1SXY",        &Tkr_1_Sxy);
    addItem("Tkr1SYY",        &Tkr_1_Syy);

    addItem("Tkr1ToTFirst",   &Tkr_1_ToTFirst);
    addItem("Tkr1ToTAve",     &Tkr_1_ToTAve);
    addItem("Tkr1ToTTrAve",   &Tkr_1_ToTTrAve);
    addItem("Tkr1ToTAsym",    &Tkr_1_ToTAsym);
    addItem("Tkr1ChisqAsym",  &Tkr_1_ChisqAsym);
    addItem("Tkr1SSDVeto",    &Tkr_1_SSDVeto);

    addItem("Tkr2Chisq",      &Tkr_2_Chisq);
    addItem("Tkr2FirstChisq", &Tkr_2_FirstChisq);

    addItem("Tkr2Hits",       &Tkr_2_Hits);
    addItem("Tkr2FirstHits",  &Tkr_2_FirstHits);
    addItem("Tkr2FirstLayer", &Tkr_2_FirstLayer);
    addItem("Tkr2LastLayer",  &Tkr_2_LastLayer);
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
    addItem("Tkr2Theta",      &Tkr_2_Theta);
    addItem("Tkr2X0",         &Tkr_2_x0);
    addItem("Tkr2Y0",         &Tkr_2_y0);
    addItem("Tkr2Z0",         &Tkr_2_z0);    

    addItem("Tkr2TkrAngle",         &Tkr_2TkrAngle); 
    addItem("Tkr2TkrHDoca",         &Tkr_2TkrHDoca); 

    // for test, uncomment these statements:
    //addItem("TkrFloat", &Tkr_float);
    //addItem("TkrInt",   &Tkr_int);

    zeroVals();

    return sc;
}

namespace {

    // coefs from Miner
    double cfThin    = 0.68;  //Set overall value by slope at 1 GeV Verticle vs 1stLayerNumber
    double cfThick   = 2.93; // Set relative value to ratio of radiation lenghts 1 : 4.33.
	                         // cfThin is close to that derived via linear regression
    double rm_hard   = 30.; 
    double rm_soft   = 130;
    double gap       = 18.; 
    double hard_frac = .7; 

    double maxToTVal =  250.;  // won't be needed after new tag of TkrDigi
    double maxPath   = 2500.;  // limit the upward propagator    
}

StatusCode TkrValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    //Tkr_float = 5.5;
    //Tkr_int   = 123;

    //offset comes from Geometry
    double z0 = m_tkrGeom->gettkrZBot();

    //special stuff here
    Tkr_1_FirstGapPlane = -1;

    double radThin  = m_tkrGeom->getAveConv(STANDARD); 
    double radThick = m_tkrGeom->getAveConv(SUPER); 
    double radTray  = m_tkrGeom->getAveRest(ALL);

    //Recover EventHeader Pointer
    //SmartDataPtr<Event::EventHeader> pEvent(m_pEventSvc, EventModel::EventHeader);

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrTrackCol>   pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol);
    if(!pTracks) return sc;

    SmartDataPtr<Event::TkrVertexCol>  pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    SmartDataPtr<Event::TkrClusterCol> pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);


    // all variable values are preset to zero. Be sure to re-initialize the ones you care about  

    double die_width = m_tkrGeom->ladderPitch();
    int nDies = m_tkrGeom->nWaferAcross();
    
    if (pTracks){   
        // Count number of tracks
        int nTracks = pTracks->size();
        Tkr_No_Tracks   = nTracks;

        if(nTracks < 1) return sc;

        // Get the first Track - it should be the "Best Track"
        Event::TkrTrackColConPtr pTrack = pTracks->begin();

        const Event::TkrTrack* track_1 = *pTrack;

        Tkr_1_Chisq        = track_1->getChiSquareSmooth();
        Tkr_1_FirstChisq   = track_1->chiSquareSegment();
        Tkr_1_FirstGaps    = track_1->getNumXFirstGaps() + track_1->getNumYFirstGaps();
        Tkr_1_Qual         = track_1->getQuality();
        Tkr_1_Type         = track_1->getStatusBits();
        Tkr_1_Hits         = track_1->getNumHits();
        Tkr_1_FirstHits    = track_1->getNumSegmentPoints();
        Tkr_1_FirstLayer   = m_tkrGeom->getLayer(track_1->front()->getTkrId());
        Tkr_1_LastLayer    = m_tkrGeom->getLayer(track_1->back()->getTkrId());
        Tkr_1_Gaps         = track_1->getNumGaps();
        Tkr_1_KalEne       = track_1->getKalEnergy(); 
        Tkr_1_ConEne       = track_1->getInitialEnergy(); 
        Tkr_1_KalThetaMS   = track_1->getKalThetaMS(); 
        Tkr_1_DifHits      = track_1->getNumXHits()-track_1->getNumYHits();

        Point  x1 = track_1->getInitialPosition();
        Vector t1 = track_1->getInitialDirection();

        Tkr_1_xdir        = t1.x();
        Tkr_1_ydir        = t1.y();
        Tkr_1_zdir        = t1.z();

        Tkr_1_x0          = x1.x();
        Tkr_1_y0          = x1.y();
        Tkr_1_z0          = x1.z();

        // theta and phi are of direction of source, hence the minus sign
        // this code replaces atan and acos used before
        Tkr_1_Phi         = (-t1).phi();
        if (Tkr_1_Phi<0.0f) Tkr_1_Phi += static_cast<float>(2*M_PI);
        Tkr_1_Theta       = (-t1).theta();

        const Event::TkrTrackParams& Tkr_1_Cov = track_1->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
        Tkr_1_Sxx         = Tkr_1_Cov.getxSlpxSlp();
        Tkr_1_Sxy         = Tkr_1_Cov.getxSlpySlp();
        Tkr_1_Syy         = Tkr_1_Cov.getySlpySlp();
        double sinPhi     = sin(Tkr_1_Phi);
        double cosPhi     = cos(Tkr_1_Phi);
        Tkr_1_ThetaErr      = t1.z()*t1.z()*sqrt(cosPhi*cosPhi*Tkr_1_Sxx + 
            2.*sinPhi*cosPhi*Tkr_1_Sxy + sinPhi*sinPhi*Tkr_1_Syy); 
        Tkr_1_PhiErr        = (-t1.z())*sqrt(sinPhi*sinPhi*Tkr_1_Sxx - 
            2.*sinPhi*cosPhi*Tkr_1_Sxy + cosPhi*cosPhi*Tkr_1_Syy);
        Tkr_1_ErrAsym     = fabs(Tkr_1_Sxy/(Tkr_1_Sxx + Tkr_1_Syy));
        Tkr_1_CovDet      = sqrt(Tkr_1_Sxx*Tkr_1_Syy-Tkr_1_Sxy*Tkr_1_Sxy)*Tkr_1_zdir*Tkr_1_zdir;

        Tkr_TrackLength = -(Tkr_1_z0-z0)/Tkr_1_zdir;

        double z_dist    = fabs((m_tkrGeom->trayHeight()+3.)/t1.z()); 
        double x_twr = globalToLocal(x1.x(), m_towerPitch, m_xNum);
        double y_twr = globalToLocal(x1.y(), m_towerPitch, m_yNum);

        double x_prj = x_twr - t1.x()*z_dist;
        double y_prj = y_twr - t1.y()*z_dist; 

        Tkr_1_TwrEdge    = (fabs(x_twr) > fabs(y_twr)) ? fabs(x_twr) : fabs(y_twr);
        Tkr_1_PrjTwrEdge = (fabs(x_prj) > fabs(y_prj)) ? fabs(x_prj) : fabs(y_prj);
        Tkr_1_TwrEdge    = m_towerPitch/2. - Tkr_1_TwrEdge;
        Tkr_1_PrjTwrEdge = m_towerPitch/2. - Tkr_1_PrjTwrEdge;

        // New section go compute gap lengths in tracker and cal
        double x_slope   = (fabs(t1.x()) > .0001)? t1.x():.00001;
        double s_x       = (sign(t1.x())*m_towerPitch/2. - x_twr)/x_slope; 
        double y_slope   = (fabs(t1.y()) > .0001)? t1.y():.00001;
        double s_y       = (sign(t1.y())*m_towerPitch/2. - y_twr)/y_slope;

        Tkr_1_TwrGap = 0.; 
        if(s_x < s_y) { // Goes out x side of CAL Module
            if(x1.z() - z0 + s_x*t1.z() > 0) {
                double s_max = fabs((x1.z()-z0)/t1.z());
                Tkr_1_TwrGap = gap/fabs(x_slope);
                if((Tkr_1_TwrGap + s_x)> s_max ) Tkr_1_TwrGap = s_max-s_x;
            }
        }
        else {          // Goes out y side
            if(x1.z() - z0 + s_y*t1.z() > 0) {
                double s_max = fabs((x1.z()-z0)/t1.z());
                Tkr_1_TwrGap = gap/fabs(y_slope);
                if((Tkr_1_TwrGap + s_y)> s_max ) Tkr_1_TwrGap = s_max-s_y;
            }
        }

        // SSD Die loaction and edge... 
        double x_die = globalToLocal(x_twr, die_width, nDies);
        double y_die = globalToLocal(y_twr, die_width, nDies);


        Tkr_1_DieEdge  = (fabs(x_die) > fabs(y_die)) ? fabs(x_die) : fabs(y_die);
        Tkr_1_DieEdge  = die_width/2. - Tkr_1_DieEdge; 

        // Section to dig out the TOT information
        double first_ToT = 0.; 
        double last_ToT  = 0.; 
        double min_ToT   = maxToTVal; 
        double max_ToT   = 0.;  
        int    hit_counter = 0; 
        double chisq_first = 0.;
        double chisq_last  = 0.; 
        Event::TkrTrackHitVecConItr pHit = track_1->begin();
        int gapId = -1;
        bool gapFound = false;
        while(pHit != track_1->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            unsigned int bits = hit->getStatusBits();
            // check if hit is in an ssd
            if ( !gapFound && !(bits & Event::TkrTrackHit::HITISSSD)) {
                Point  gapPos = hit->getPoint(Event::TkrTrackHit::PREDICTED);
                Tkr_1_GapX = gapPos.x();
                Tkr_1_GapY = gapPos.y();
                //TkrId is good!
                gapId = m_tkrGeom->getPlane(hit->getTkrId());
                gapFound = true;
            }
            //plane--;
            if (!(bits & Event::TkrTrackHit::HITONFIT)) continue;
            const Event::TkrCluster* cluster = hit->getClusterPtr();
            int size =  (int) (const_cast<Event::TkrCluster*>(cluster))->size();
            // get the local slopes
            double slope  = fabs(hit->getMeasuredSlope(Event::TkrTrackHit::SMOOTHED));
            double slope1 = fabs(hit->getNonMeasuredSlope(Event::TkrTrackHit::SMOOTHED));

            // theta1 is the projected angle across the strip
            double theta1       = atan(slope1);

            double aspectRatio = 0.228/0.400;
            double totMax      =  250.;   // counts
            double threshold   =  0.25;   // Mips
            double countThreshold = 15; // counts
            double normFactor  =  1./53.;

            double mips = cluster->getMips();

            double tot = cluster->ToT();
            if(tot>=totMax) tot = totMax;
            double path1 = 1.0;

            // get the path length for the hit
            // tries to get the average
            // the calculation is part analytic, part approximation and part fudge.
            //   more work is definitely in order!

            // theta1 first
            if (tot>=totMax) { tot = normFactor*(totMax+countThreshold); }
            else {
                double costh1 = cos(theta1);
                if (size==1) {
                    double sinth1 = sin(theta1);
                    if (slope1< aspectRatio) {
                        path1 = (1./costh1*(aspectRatio-slope1) + 
                            (1/costh1 - 0.5*threshold)*(2*threshold*sinth1))
                            /(aspectRatio - slope1 + 2*threshold*sinth1);
                    } else if (slope1<aspectRatio/(1-2.*threshold*costh1)) {
                        path1 = 1; //1/costh1 - threshold*costh1;
                    } else { 
                        path1 = 1;
                    }
                }
                else if (size==2) {
                    if (slope1<aspectRatio/(1.-threshold*costh1)) {
                        path1 = 0.75/costh1 -0.5*threshold;
                    } else if (slope1<2.*aspectRatio/(1.-2*threshold*costh1)) { 
                        path1 = aspectRatio/sin(theta1);
                    } else {
                        path1 = 1.0;
                    }
                } else {
                    if(slope1>aspectRatio/(1.- 2.*threshold*costh1)) {
                        path1 = aspectRatio/sin(theta1);
                    } else {
                        path1 = 1.0;
                    }
                }
                double factor = path1*costh1*slope;
                double path2 = sqrt(path1*path1 + factor*factor);
                mips /= path2;
            }

            if(mips > max_ToT) max_ToT = mips; 
            if(mips < min_ToT) min_ToT = mips; 
            hit_counter++;  
            if (hit_counter==1) Tkr_1_ToTFirst = mips;
            Tkr_1_ToTAve += mips;
            if(hit_counter < 3) {
                first_ToT += mips;
                chisq_first += hit->getChiSquareSmooth();
            }
            if(hit_counter > Tkr_1_Hits - 2){
                last_ToT += mips;
                chisq_last += hit->getChiSquareSmooth();
            }
        }
        Tkr_1_ToTTrAve = (Tkr_1_ToTAve - max_ToT - min_ToT)/(Tkr_1_Hits-2.);
        Tkr_1_ToTAve /= Tkr_1_Hits;
        Tkr_1_ToTAsym = (last_ToT - first_ToT)/(first_ToT + last_ToT);
        Tkr_1_FirstGapPlane = gapId; 


        // Chisq Asymmetry - Front vs Back ends of tracks
        Tkr_1_ChisqAsym = (chisq_last - chisq_first)/(chisq_last + chisq_first);

        m_G4PropTool->setStepStart(x1, -t1); //Note minus sign - swim backwards towards ACD

        int topPlane = m_tkrGeom->numPlanes()-1; 
        double topOfTkr = m_tkrGeom->getPlaneZ(topPlane) + 1.0;
        double arc_min = fabs((topOfTkr-x1.z())/t1.z());
        arc_min = std::min( arc_min, maxPath); 
        m_G4PropTool->step(arc_min);  
        int numSteps = m_G4PropTool->getNumberSteps();

        idents::VolumeIdentifier volId;
        idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();

        for(int istep = 1; istep < numSteps; ++istep) { // Note: skip the first vol - as this is the head SSD
            volId = m_G4PropTool->getStepVolumeId(istep);
            volId.prepend(prefix);
            Point x_step       = m_G4PropTool->getStepPosition(istep); 
            if((x_step.z()-x1.z()) < 10.0) continue; 
            if(x_step.z() > topOfTkr || !m_tkrGeom->isInActiveLAT(x_step) ) break; 

            // check that it's really a TKR hit (probably overkill)
            if(volId.size() != 9) continue; 
            if(!(volId[0]==0 && volId[3]==1)) continue; // !(LAT && TKR)
            if(volId[6]> 1) continue;  //It's a converter!  

            Tkr_1_SSDVeto += 1.; 
        }

        if(nTracks > 1) {
            pTrack++;
            const Event::TkrTrack* track_2 = *pTrack;

            Tkr_2_Chisq        = track_2->getChiSquareSmooth();
            Tkr_2_FirstChisq     = track_2->chiSquareSegment();
            Tkr_2_FirstGaps      = track_2->getNumXFirstGaps() + track_2->getNumYFirstGaps();
            Tkr_2_Qual         = track_2->getQuality();
            Tkr_2_Type         = track_2->getStatusBits();
            Tkr_2_Hits         = track_2->getNumHits();
            Tkr_2_FirstHits    = track_2->getNumSegmentPoints();
            Tkr_2_FirstLayer   = m_tkrGeom->getLayer(track_2->front()->getTkrId());
            Tkr_2_LastLayer    = m_tkrGeom->getLayer(track_2->back()->getTkrId());
            Tkr_2_Gaps         = track_2->getNumGaps();
            Tkr_2_KalEne       = track_2->getKalEnergy(); 
            Tkr_2_ConEne       = track_2->getInitialEnergy(); 
            Tkr_2_KalThetaMS   = track_2->getKalThetaMS(); 
            Tkr_2_DifHits      = track_2->getNumXHits()-track_2->getNumYHits();

            Point  x2 = track_2->getInitialPosition();
            Vector t2 = track_2->getInitialDirection();
            Tkr_2_xdir       = t2.x();
            Tkr_2_ydir       = t2.y();
            Tkr_2_zdir       = t2.z();

            // this replaces atan used before
            Tkr_2_Phi         = (-t2).phi();
            if (Tkr_2_Phi<0.0) Tkr_2_Phi += static_cast<float>(2*M_PI);
            Tkr_2_Theta       = (-t2).theta();

            Tkr_2_x0         = x2.x();
            Tkr_2_y0         = x2.y();
            Tkr_2_z0         = x2.z();

            double x_twr = globalToLocal(x2.x(), m_towerPitch, m_xNum);
            double y_twr = globalToLocal(x2.y(), m_towerPitch, m_yNum);
            double x_prj = x_twr - t2.x()*z_dist;
            double y_prj = y_twr - t2.y()*z_dist; 

            Tkr_2_TwrEdge    = (fabs(x_twr) > fabs(y_twr)) ? fabs(x_twr) : fabs(y_twr);
            Tkr_2_PrjTwrEdge = (fabs(x_prj) > fabs(y_prj)) ? fabs(x_prj) : fabs(y_prj);
            Tkr_2_TwrEdge    = m_towerPitch/2. - Tkr_2_TwrEdge;
            Tkr_2_PrjTwrEdge = m_towerPitch/2. - Tkr_2_PrjTwrEdge;

            double x_die = globalToLocal(x_twr, die_width, nDies);
            double y_die = globalToLocal(y_twr, die_width, nDies);

            Tkr_2_DieEdge  = (fabs(x_die) > fabs(y_die)) ? fabs(x_die) : fabs(y_die);
            Tkr_2_DieEdge  = die_width/2. - Tkr_2_DieEdge; 

            Tkr_2TkrAngle = acos(t1*t2);  
            Point x2p  = x2 + ((x1.z()-x2.z())/t2.z())*t2;
            Point x20  = x2 - (x2.z()/t2.z())*t2;
            Point x10  = x1 - (x1.z()/t1.z())*t1;
            double doca_plane = (x2p-x1).mag();
            double doca_0     = (x20-x10).mag();
            if(doca_plane > doca_0) Tkr_2TkrAngle *= -1.; 
            Tkr_2TkrHDoca = -doca_plane*t1.z();
        }

        Tkr_Sum_KalEne    = Tkr_1_KalEne+Tkr_2_KalEne; 
        Tkr_Sum_ConEne    = Tkr_1_ConEne+Tkr_2_ConEne;      


        // Computation of the tracker contribution to the total energy 
        double costh = fabs(t1.z());
        double secth = 1./costh;
        arc_min = (x1.z() - m_tkrGeom->calZTop())*secth; 
        m_G4PropTool->setStepStart(x1, t1);
		m_G4PropTool->step(arc_min);
		double z_present = x1.z();

        // Compute the sum-of radiation_lengths x Hits in each layer
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

        int firstPlane = m_tkrGeom->getPlane(track_1->front()->getTkrId()); 
        int firstLayer = m_tkrGeom->getLayer(firstPlane); 


        // doesn't work for reverse-found tracks
        for(int ilayer = firstLayer; ilayer>=0; --ilayer) {
            double xms = 0.;
            double yms = 0.;

            if(ilayer <firstLayer) {
                HepMatrix Q = m_G4PropTool->getMscatCov(arc_len, Tkr_Sum_ConEne/2.);
                xms = Q(1,1);
                yms = Q(3,3);
                radlen = m_G4PropTool->getRadLength(arc_len); 
            }
            double xSprd = 80.*secth; 
            double ySprd = 80.*secth;  
            double halfTray = 0.5*m_tkrGeom->trayWidth();
            xSprd = std::min(xSprd, halfTray);
            ySprd = std::min(ySprd, halfTray);

            // Assume location of shower center is given by 1st track
            Point x_hit = x1 + arc_len*t1;
            int numHits = pQueryClusters->numberOfHitsNear(ilayer, xSprd, ySprd, x_hit, t1);

			//  Look for extra hits around track head
            if(ilayer == firstLayer) {
                double xRgn = 30.*secth;
                double yRgn = 30.*secth;
                Tkr_HDCount = pQueryClusters->numberOfUUHitsNear(ilayer, xRgn, yRgn, x_hit);
            }

            double layer_edge = towerEdge(x_hit);

            double delta_rad= radlen-radlen_old;
            double thisRad;
            //A bit cleaner
            switch (m_tkrGeom->getLayerType(ilayer)) {
            case STANDARD:
                thisRad = radThin;
                thin_hits += numHits;
                break;
            case SUPER:
                thisRad = radThick;
                thick_hits += numHits;
            case NOCONV:
                thisRad = 0.0;
                blank_hits += numHits;
                break;
            default:
                break;
            }
            if (ilayer==firstLayer) {
                // on first layer, add in 1/2 if the 1st plane is the top of a layer
                if(m_tkrGeom->isTopPlaneInLayer(firstPlane)) {delta_rad = 0.5*thisRad*secth;}
            } else {
                // on subseqent layers, make sure that there is a minimum radiator
                if(delta_rad*costh < thisRad) {
                    delta_rad = (radTray + thisRad)*secth;
                }
            }

            total_hits       += numHits; 
            ave_edge         += layer_edge*delta_rad; 
            rad_len_sum      += delta_rad;

            // Increment arc-length
            if(ilayer==0) break;

			double z_next = m_tkrGeom->getLayerZ(ilayer-1);
            double deltaZ = z_present - z_next;
			z_present = z_next;

            arc_len += fabs( deltaZ/t1.z()); 
            radlen_old = radlen; 
        }
 
        Tkr_Energy     = (cfThin*thin_hits + cfThick*thick_hits)/std::max(costh, .2);

		// The following flattens the cos(theta) dependence.  Anomolous leakage for widely spaced
		// samples?  
        Tkr_Energy_Corr= Tkr_Energy*(1.+ .0012*(Tkr_1_FirstLayer-1)*(Tkr_1_FirstLayer-1)) 
			                        *(1 + .3*std::max((4-Tkr_1_FirstLayer),0.f));

        Tkr_Total_Hits = total_hits;
        Tkr_Thin_Hits  = thin_hits;
        Tkr_Thick_Hits = thick_hits;
        Tkr_Blank_Hits = blank_hits; 
        Tkr_TwrEdge    = ave_edge/rad_len_sum; 
        Tkr_RadLength  = rad_len_sum;
    }
    return sc;
}

double TkrValsTool::towerEdge(Point pos) const
{
    double edge = 0.; 
    double x = pos.x();
    double y = pos.y();

    double x_twr = globalToLocal(x, m_towerPitch, m_xNum);
    double y_twr = globalToLocal(y, m_towerPitch, m_yNum);

    edge = 0.5*m_towerPitch - std::max(fabs(x_twr),fabs(y_twr));
    return edge;
}