/** @file TkrValsTool.cxx
@brief Calculates the Tkr analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrValsTool.cxx,v 1.41 2004/09/10 04:59:02 lsrea Exp $
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
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitPlane.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrFitPar.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrFitMatrix.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
//#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GaudiKernel/IToolSvc.h"

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
    ITkrGeometrySvc*       pTkrGeoSvc;
    /// GlastDetSvc used for access to detector info
    IGlastDetSvc*         m_detSvc; 
    /// pointer to queryclusterstool
    ITkrQueryClustersTool* pQueryClusters;
    /// 
    IPropagatorSvc* m_propSvc;

    IKalmanParticle*       pKalParticle;

    IPropagator*           m_G4PropTool;    

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
	double Tkr_1_FirstGapLayer; 
	double Tkr_1_GapX;
	double Tkr_1_GapY;
    double Tkr_1_FirstGaps; 
    double Tkr_1_Hits;
    double Tkr_1_FirstHits;
    double Tkr_1_FirstLayer; 
    double Tkr_1_LastLayer; 

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
    double Tkr_1_ThetaErr;
    double Tkr_1_PhiErr;
    double Tkr_1_ErrAsym;
    double Tkr_1_CovDet;
    double Tkr_1_ToTFirst;
    double Tkr_1_ToTAve;
    double Tkr_1_ToTTrAve;
    double Tkr_1_ToTAsym;
    double Tkr_1_ChisqAsym;
    double Tkr_1_SSDVeto; 

    //Second Track Specifics
    double Tkr_2_Chisq;
    double Tkr_2_FirstChisq;
    double Tkr_2_FirstGaps; 
    double Tkr_2_Qual;
    double Tkr_2_Type;
    double Tkr_2_Hits;
    double Tkr_2_FirstHits;
    double Tkr_2_FirstLayer; 
    double Tkr_2_LastLayer; 

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
    double Tkr_2_Theta;
    double Tkr_2_x0;
    double Tkr_2_y0;
    double Tkr_2_z0;

    double Tkr_2TkrAngle;
    double Tkr_2TkrHDoca;

    // here's some test stuff... if it works for a couple it will work for all
    float Tkr_float;
    int   Tkr_int;
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

        m_towerPitch = pTkrGeoSvc->towerPitch();
        m_xNum       = pTkrGeoSvc->numXTowers();
        m_yNum       = pTkrGeoSvc->numYTowers();

        // find GlastDevSvc service
        if (service("GlastDetSvc", m_detSvc, true).isFailure()){
            log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
            return StatusCode::FAILURE;
        }

        // pick up the chosen propagator
        if (service("GlastPropagatorSvc", m_propSvc, true ).isFailure()) {
            log << MSG::ERROR << "Couldn't find the GlastPropagatorSvc!" << endreq;
            return fail;
        }  
        pKalParticle = m_propSvc->getPropagator();

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
    addItem("Tkr1LastLayer",  &Tkr_1_LastLayer);
    addItem("Tkr1DifHits",    &Tkr_1_DifHits);

    addItem("Tkr1Gaps",       &Tkr_1_Gaps);
    addItem("Tkr1FirstGapLayer",&Tkr_1_FirstGapLayer);
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

    double maxPath = 2500.; // limit the upward propagator    
}

StatusCode TkrValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    Tkr_float = 5.5;
    Tkr_int   = 123;

    //placeholder for offset
    double z0 = 0.0;

    double radThin  = pTkrGeoSvc->getAveConv(STANDARD); // was 0.03
    double radThick = pTkrGeoSvc->getAveConv(SUPER);    // was 0.18
    double radTray  = pTkrGeoSvc->getAveRest(ALL);      // was 0.015

    ITkrToTSvc* pToTSvc = pTkrGeoSvc->getTkrToTSvc();

    //Recover EventHeader Pointer
    //SmartDataPtr<Event::EventHeader> pEvent(m_pEventSvc, EventModel::EventHeader);

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    //SmartDataPtr<Event::TkrVertexCol>     pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    SmartDataPtr<Event::TkrClusterCol> pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    if(!pTracks) return StatusCode::FAILURE;

    // all variable values are preset to zero. Be sure to re-initialize the ones you care about  
    
    //Make sure we have valid reconstructed data

    int nNoConv = pTkrGeoSvc->numNoConverter();
    int nThick  = pTkrGeoSvc->numSuperGlast();
    int nThin   = pTkrGeoSvc->numLayers() - nThick - nNoConv;

    double die_width = pTkrGeoSvc->ladderNStrips()*pTkrGeoSvc->siStripPitch() +
        + 2.*pTkrGeoSvc->siDeadDistance() + pTkrGeoSvc->ladderGap();
    int nDies = pTkrGeoSvc->nWaferAcross();

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
        Tkr_1_LastLayer    = track_1->getLayer(Event::TkrFitTrackBase::End);
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


        // theta and phi are of direction of source, hence the minus sign
        // this code replaces atan and acos used before
        Tkr_1_Phi         = (-t1).phi();
        if (Tkr_1_Phi<0.0) Tkr_1_Phi += 2*M_PI;
        Tkr_1_Theta       = (-t1).theta();

        Event::TkrFitMatrix  Tkr_1_Cov = track_1->getTrackCov();
        Tkr_1_Sxx         = Tkr_1_Cov.getcovSxSx();
        Tkr_1_Sxy         = Tkr_1_Cov.getcovSxSy();
        Tkr_1_Syy         = Tkr_1_Cov.getcovSySy();
        double sinPhi     = sin(Tkr_1_Phi);
        double cosPhi     = cos(Tkr_1_Phi);
        Tkr_1_ThetaErr      = t1.z()*t1.z()*sqrt(cosPhi*cosPhi*Tkr_1_Sxx + 
            2.*sinPhi*cosPhi*Tkr_1_Sxy + sinPhi*sinPhi*Tkr_1_Syy); 
        Tkr_1_PhiErr        = (-t1.z())*sqrt(sinPhi*sinPhi*Tkr_1_Sxx - 
            2.*sinPhi*cosPhi*Tkr_1_Sxy + cosPhi*cosPhi*Tkr_1_Syy);
        Tkr_1_ErrAsym     = fabs(Tkr_1_Sxy/(Tkr_1_Sxx + Tkr_1_Syy));
        Tkr_1_CovDet      = sqrt(Tkr_1_Sxx*Tkr_1_Syy-Tkr_1_Sxy*Tkr_1_Sxy)*Tkr_1_zdir*Tkr_1_zdir;

        Tkr_TrackLength = -(Tkr_1_z0-z0)/Tkr_1_zdir;

        double z_dist    = fabs((pTkrGeoSvc->trayHeight()+3.)/t1.z()); 
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
        double totMax      =  pToTSvc->getMaxToT();   // counts
        double first_ToT = 0.; 
        double last_ToT  = 0.; 
        double min_ToT   = totMax; 
        double max_ToT   = 0.;  
        int    hit_counter = 0; 
        double chisq_first = 0.;
        double chisq_last  = 0.; 
        Event::TkrFitPlaneConPtr pln_pointer = track_1->begin();
        //for the ToT path correction; use the directions of each hit
        double slopeX;
        double slopeY;

		int gapId = -1; 
		int lastPlane = -1; 
        while(pln_pointer != track_1->end()) {
            Event::TkrFitPlane plane = *pln_pointer;

			int thisLayer = plane.getIDPlane();
			int thisView = Event::TkrCluster::viewToInt(plane.getProjection());
            int tray, botTop;
			pTkrGeoSvc->layerToTray(thisLayer, thisView, tray, botTop);
            int thisPlane = pTkrGeoSvc->trayToPlane(tray, botTop);
            // first hit
			if(lastPlane < 0) { lastPlane = thisPlane;
			} else {
				if(gapId < 0 && lastPlane+1 != thisPlane){
					gapId = lastPlane+1;
					Event::TkrFitPlane lastFitPlane = *(--pln_pointer);
					pln_pointer++;
					Point lastPoint = lastFitPlane.getPoint(Event::TkrFitHit::FIT);
					Event::TkrFitHit lastHit = lastFitPlane.getHit(Event::TkrFitHit::FIT);
					double xSlope = lastHit.getPar().getXSlope();
					double ySlope = lastHit.getPar().getYSlope();
					Vector localDir = Vector(-xSlope,-ySlope,-1.).unit();
					Ray localSeg(lastPoint, localDir);
					int view; 
					int layer; 
					pTkrGeoSvc->planeToLayer (gapId, layer, view);
					double gapZ = pTkrGeoSvc->getReconLayerZ(layer,view);
					double arcLen = (gapZ-lastPoint.z())/localDir.z();
					Point gapPoint = localSeg.position(arcLen);
					Tkr_1_GapX = gapPoint.x();
					Tkr_1_GapY = gapPoint.y();
				}
			}
			lastPlane = thisPlane;

            int hit_Id = plane.getIDHit();
            Event::TkrCluster* cluster = pClusters->getHit(hit_Id);
            Event::TkrCluster::view v = cluster->v();
            int size =  (int) cluster->size();

            // get the local slopes
            Event::TkrFitPar par = (plane.getHit(Event::TkrFitHit::SMOOTH )).getPar();
            slopeX = fabs(par.getXSlope());
            slopeY = fabs(par.getYSlope());
            double slope        = (v==Event::TkrCluster::X) ? slopeY : slopeX;
            double slope1       = (v==Event::TkrCluster::X) ? slopeX : slopeY;

            // theta is the projected angle along the strip
            // theta1 is the projected angle across the strip
            double theta1       = atan(slope1);

            double aspectRatio = 0.228/0.400;
            double threshold   =  0.25;   // Mips (for old-style)
            double countThreshold = 15; // counts (for old-style)
            double normFactor  =  1./53.; // (for old-style)

            double mips;
            double rawToT = cluster->getRawToT();
            if (cluster->getVersion()==0) {
                mips = normFactor*(std::min(totMax,rawToT)+countThreshold);
            } else {
                mips = cluster->getMips();
            }

            double path1 = 1.0;

            // get the path length for the hit
            // the calculation is part analytic, part approximation and part fudge.
            //   more work is definitely in order!

            // theta1 first
            if(rawToT<totMax) {
                double costh1 = cos(theta1);
                if (size==1) {
                    if (slope1< aspectRatio) {
                        path1 = (1./costh1*(aspectRatio-slope1) + 
                            (1/costh1 - 0.5*threshold)*(2*threshold*sin(theta1)))
                            /(aspectRatio - slope1 + 2*threshold*sin(theta1));
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
                chisq_first += plane.getDeltaChiSq(Event::TkrFitHit::SMOOTH);
            }
            if(hit_counter > Tkr_1_Hits - 2){
                last_ToT += mips;
                chisq_last += plane.getDeltaChiSq(Event::TkrFitHit::SMOOTH);
            }
            pln_pointer++;

        }
        Tkr_1_ToTTrAve = (Tkr_1_ToTAve - max_ToT - min_ToT)/(Tkr_1_Hits-2.);
        Tkr_1_ToTAve /= Tkr_1_Hits;
        Tkr_1_ToTAsym = (last_ToT - first_ToT)/(first_ToT + last_ToT);
		Tkr_1_FirstGapLayer = gapId; 

        // Chisq Asymetry - Front vs Back ends of tracks
        Tkr_1_ChisqAsym = (chisq_last - chisq_first)/(chisq_last + chisq_first);

        m_G4PropTool->setStepStart(x1, -t1); //Note minus sign - swim backwards towards ACD

        double topOfTkr = pTkrGeoSvc->getReconLayerZ(0) + 2.0; // higher than first silicon layer
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
            if(x_step.z() > topOfTkr || !pTkrGeoSvc->isInActiveLAT(x_step) ) break; 
            //if(x_step.z() > topOfTkr || fabs(x_step.x()) > xEdge || fabs(x_step.y()) > yEdge) break; 

            // check that it's really a TKR hit (probably overkill)
            if(volId.size() != 9) continue; 
            if(!(volId[0]==0 && volId[3]==1)) continue; // !(LAT && TKR)
            if(volId[6]> 1) continue;  //It's a converter!  

            Tkr_1_SSDVeto += 1.; 
        }

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
            Tkr_2_FirstHits    = track_2->getNumSegmentPoints();
            Tkr_2_FirstLayer   = track_2->getLayer();
            Tkr_2_LastLayer    = track_2->getLayer(Event::TkrFitTrackBase::End);
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

            // this replaces atan used before
            Tkr_2_Phi         = (-t2).phi();
            if (Tkr_2_Phi<0.0) Tkr_2_Phi += 2*M_PI;
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
        arc_min = (x1.z() - pTkrGeoSvc->calZTop())/costh; 
        pKalParticle->setStepStart(x1, t1, arc_min);

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

            //int outside = 0; 
            //int iView   = 0;
            double layer_edge = towerEdge(x_hit);

            double in_frac_soft = containedFraction(x_hit, gap, rm_soft, costh, Tkr_1_Phi);
            if(in_frac_soft < .01) in_frac_soft = .01; 

            double in_frac_hard = containedFraction(x_hit, gap, rm_hard, costh, Tkr_1_Phi);
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

            double ene_corr = corr_factor*delta_rad; //*ene; 
            tracker_ene_corr += ene_corr;
            total_hits       += numHits; 
            ave_edge         += layer_edge*delta_rad; 
            rad_len_sum      += delta_rad;

            // Increment arc-length
            int nextPlane = iplane+1;
            if (iplane==max_planes-1) nextPlane--;
            double deltaZ = pTkrGeoSvc->getReconLayerZ(iplane) -
                pTkrGeoSvc->getReconLayerZ(nextPlane);
            arc_len += fabs( deltaZ/t1.z()); 
            radlen_old = radlen; 
        }
        // Coef's from Linear Regression analysis in Miner
        // defined in anonymous namespace above
		// Note: there is a corner of phase space where this can get much too large
		//  (a check would be to compare with TkrKalEneSum - if its larger, pick
		//   TkrKalEneSum but maybe not -  this will happen in mis-tracked events)
        Tkr_Energy     = (cfThin*thin_hits + cfThick*thick_hits + cfNoConv*blank_hits
			             + cfRadLen*std::min(rad_len_sum, 3.0)
			             + cfZ*(x1.z()-z0))/std::max(costh, .4);
 
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

   // std::cout << "Tkr_KalEne: " << Tkr_1_KalEne << " " << Tkr_2_KalEne << std::endl 
   //              <<" Tkr_KalThetaMS " <<  Tkr_1_KalThetaMS << " " << Tkr_2_KalThetaMS << std::endl;

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

double  TkrValsTool::containedFraction(Point pos, double gap,  
                                       double r, double costh, double phi) const
{
    // Get the projected angles for the gap
    double tanth = sqrt(1.-costh*costh)/costh;
    double gap_x = gap - 35.*sin(phi)*tanth;
    double x = pos.x();
    double y = pos.y();
    if(gap_x < 5.) gap_x = 5.;
    double gap_y = gap - 35.*cos(phi)*tanth;
    if(gap_y < 5.) gap_y = 5.; 

    // X Edges
    double x_twr = globalToLocal(x, m_towerPitch, m_xNum);
    double edge = m_towerPitch/2. - fabs(x_twr);
    double r_frac_plus = (edge-gap_x/2.)/r; 
    double angle_factor = sin(phi)*(1./costh - 1.);
    double in_frac_x  =  circleFractionSimpson(r_frac_plus, angle_factor);
    if (x>pTkrGeoSvc->getLATLimit(0,LOW)+0.5*m_towerPitch
        && x<pTkrGeoSvc->getLATLimit(0,HIGH)-0.5*m_towerPitch) {
        double r_frac_minus = (edge + gap_x/2.)/r;
        in_frac_x += circleFractionSimpson(-r_frac_minus, angle_factor);
    }

    // Y Edges
    double y_twr = globalToLocal(y, m_towerPitch, m_yNum);
    edge = m_towerPitch/2. - fabs(y_twr);
    r_frac_plus = (edge-gap_y/2.)/r; 
    angle_factor = cos(phi)*(1./costh - 1.);
    double in_frac_y  =  circleFractionSimpson(r_frac_plus, angle_factor);
    if (y>pTkrGeoSvc->getLATLimit(1,LOW)+0.5*m_towerPitch
        && y<pTkrGeoSvc->getLATLimit(1,HIGH)-0.5*m_towerPitch) {
        double r_frac_minus = (edge + gap_y/2.)/r;
        in_frac_y += circleFractionSimpson(-r_frac_minus, angle_factor);
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
