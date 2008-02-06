/** @file TkrValsTool.cxx
@brief Calculates the Tkr analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrValsTool.cxx,v 1.82.24.1 2008/02/06 20:47:11 heather Exp $
*/
//#define PRE_CALMOD 1

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
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "geometry/Ray.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrFlagHitsTool.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GaudiKernel/IToolSvc.h"
#include "Doca.h"

#include <cstring>


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
    double m_activeWidth;
    bool   m_useNew;
    bool   m_enableVetoDiagnostics;

    // some pointers to services

    /// pointer to tracker geometry
    ITkrGeometrySvc*       m_tkrGeom;
    /// GlastDetSvc used for access to detector info
    IGlastDetSvc*         m_detSvc; 
    /// pointer to queryclusterstool
    ITkrQueryClustersTool* pQueryClusters;
    /// pointer to flagHitsTool
    ITkrFlagHitsTool* pFlagHits;
    /// 
    IPropagatorSvc* m_propSvc;

    IPropagator*           m_G4PropTool; 

    // properties
    double m_minVetoError;
    double m_maxVetoError;
    double m_vetoNSigma;


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
    float Tkr_SurplusHitRatio;
    float Tkr_SurplusHCInside;
    float Tkr_UpstreamHC;
    float TkrDispersion;

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
    float Tkr_1_SSDVetoOld; 
    float Tkr_1_SSDVeto; 

    // for SSDVeto Diagnostics
    int   Tkr_1_VetoTrials;
    int   Tkr_1_VetoHitFound;
    int   Tkr_1_VetoUnknown;
    int   Tkr_1_VetoPlaneCrossed;
    int   Tkr_1_VetoTower;
    int   Tkr_1_VetoGapCorner;
    //double Tkr_1_MinGapDistance;
    //double Tkr_1_MaxGapDistance;
    int   Tkr_1_VetoGapEdge;
    int   Tkr_1_VetoBadCluster;
    int   Tkr_1_VetoDeadPlane;
    int   Tkr_1_VetoTruncated;

    float Tkr_1_CoreHC;
    float Tkr_1_LATEdge;

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

    /*
    float Tkr_float;
    int   Tkr_int;
    float    Tkr_1Pos[3];
    double   Tkr_doubles[2];
    int      Tkr_ints[5];
    unsigned int Tkr_uInts[7];
    char Tkr_string[20];
    */
};

namespace 
{
    double interpolate ( double xi, double xmin, double xmax, double ymin, double ymax) 
    {
        double step = (xi-xmin)/(xmax-xmin);
        return ymin + (ymax-ymin)*step;
    }

    int xPosIdx = Event::TkrTrackParams::xPosIdx;
    int yPosIdx = Event::TkrTrackParams::yPosIdx;
    int xSlpIdx = Event::TkrTrackParams::xSlpIdx;
    int ySlpIdx = Event::TkrTrackParams::ySlpIdx;
}

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

    declareProperty("useNew", m_useNew=true);
    declareProperty("enableVetoDiagnostics", m_enableVetoDiagnostics=true);
    declareProperty("minVetoError", m_minVetoError=1.0);
    declareProperty("maxVetoError", m_maxVetoError=100000.0);
    declareProperty("vetoNSigma", m_vetoNSigma=2.0);
}

/** @page anatup_vars 
    @section tkrvalstool TkrValsTool Variables

Notes: 
- Variables called Tkr1Xxx refer to the "best" track; 
     those called Tkr2Xxx refer to the second track.
- For variables listed as Tkr[1/2]Xxx there are two versions in the ntuple, 
     one for the best and one for the second track. 
- The labels are not entierly consistent, but it's probably 
     too disruptive to fix them at this point.
     For example: TkrRadLength, TkrTrackLength, TkrTwrEdge refer to track 1. 
     Also, Tkr2Angle and Tkr2HDoca are quantities that depend on both tracks.
- The variables associated with the second track are undefined 
     if there is only one track! 
     Check TkrNumTracks before using these variables! 
     In fact check TkrNumTracks before using first-track variables, for the same reason.
- A new section of (optional) ssd-veto diagnostic variables has been added. They are
     included by default.

@subsection general General variables
     <table>
<tr><th> Variable <th> Type  <th> Description				                 
<tr><td> TkrNumTracks 	
<td>F<td>   Number of tracks found (Maximum is set by TkrRecon, currently 10) 
<tr><td> TkrSumKalEne 
<td>F<td>   Sum of Kalman energies (see TkrNKalEne, below) 
            for the two best tracks 
<tr><td> TkrSumConEne 	
<td>F<td>   Sum of the energies for the two best tracks, 
            as assigned by the patrec energy tool 
<tr><td> TkrEnergy 
<td>F<td>   Energy in tracker, as determined from linear regression analysis 
            of number of clusters 
<tr><td> TkrEnergySum 	
<td>F<td>   Deprecated 
<tr><td> TkrEnergyCorr 
<td>F<td>   TkrEnergy corrected by TkrEdgeCorr 
<tr><td> TkrEdgeCorr 	
<td>F<td>   Tracker edge correction. This may go away; it's an intermediate quantity 
<tr><td> TkrHDCount 
<td>F<td>   Number of unused clusters in top x-y layer of the best track 
            within a radius of 30 mm, corrected for track angle
            (Used in PSF analysis and background rejection) 
<tr><td> TkrTotalHits 
<td>F<td>   Deprecated. Use TkrSurplusHCInside instead
<tr><td> TkrSurplusHitsInside 
<td>F<td>   Number of clusters inside an energy- and angle-dependent cone 
            centered on the reconstructed axis of the best track and
            starting at the head of track 1. Only hits in layers with at
            least one x and one y cluster in the tower are counted.
<tr><td> TkrSurplusHitRatio
<td>F<td>   Ratio of the number of hits outside the cone to the number
            inside. See TkrSurplusHitsInside
<tr><td> TkrThinHits
<td>F<td>   Number of clusters in the above cone in the thin-converter layers 
<tr><td> TkrThickHits 
<td>F<td>   Number of clusters in the above cone in the thick-converter layers 
<tr><td> TkrBlankHits 
<td>F<td>   Number of clusters in the above cone in the no-converter layers 
<tr><td> Tkr2TkrAngle 
<td>F<td> Angle between first and second reconstructed tracks 
<tr><td> Tkr2TkrHDoca  
<td>F<td>   Distance between first and second track in the plane of the 
            first hit on the first track. 
            This is most useful if the two tracks are almost parallel, 
            in which case the usual DOCA is poorly measured.  
</table>

@subsection both Variables that exist for both best and second tracks

<table>
<tr><th> Variable <th> Type  <th> Description				                 
<tr><td> Tkr[1/2]Chisq 
<td>F<td>   Track chisquared 
<tr><td> Tkr[1/2]FirstChisq  
<td>F<td>   Track chisquared for first Tkr[1/2]FirstHits layers  
<tr><td> Tkr[1/2]Hits  
<td>F<td>   Number of clusters in track  
<tr><td> Tkr[1/2]FirstHits  
<td>F<td>   Number of initial track hits used to determine the starting direction  
<tr><td> Tkr[1/2][First/Last]Layer  
<td>F<td>   [First/Last] layer in track  (layer 0 is the bottom of the tracker)
<tr><td> Tkr[1/2]DifHits  
<td>F<td>   Difference between the number of x and y clusters associated with track  
<tr><td> Tkr[1/2]Gaps  
<td>F<td>   Total number of gaps in track  
<tr><td> Tkr1FirstGapPlane  
<td>F<td>   plane number of first gap on track 1  
             (This and the following X,Y pair can be used to find dead strips)
<tr><td> Tkr1[X/Y]Gap  
<td>F<td>   [x/y] location of first gap on track 1  
<tr><td> Tkr[1/2]FirstGaps  
<td>F<td>   Number of gaps in first Tkr1FirstHits layers on track  
<tr><td> Tkr[1/2]Qual  
<td>F<td>   Track "quality": depends on the number of clusters and chisquared of the track. 
             Maximum is currently 64, can be negative if chisqared gets large. 
             This is used primarily to order the tracks during patrec. 
             <strong>It's not a good idea to cut on this variable!</strong>  
<tr><td> Tkr[1/2]Type  
<td>F<td>   These are the status bits from the trackign, containing information
             about how the track was found and fitted.   
             See TkrTrack.h in the Event package for the current description.     
             As of GlastRelease v7r2, the status word bits organized as follows:

@verbatim

       |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
        [ Pat Rec Info  ]  [Pass ][ E-Loss] [ Track Energy ]  [Track Fit Status]

        FOUND    = 0x0001,  //Set if track has been "found" by pat rec
        FILTERED = 0x0002,  //Set if track fit filter stage has been run
        SMOOTHED = 0x0004,  //Set if track fit smoother has been run
        REVFILTR = 0x0008,  //Set if track has been reverse-filtered
        CALENERGY= 0x0010,  //Set if track energy from raw calorimeter info
        LATENERGY= 0x0020,  //Set if track energy from TKR+CAL constrained
        USERENERGY= 0x0040, //Set if track energy set by user
        MCENERGY = 0x0080,  //Set if energy from users or from MC truth
        RADELOSS = 0x0100,  //Set if radiative energy loss used (e+/e- fitting)
        MIPELOSS = 0x0200,  //Set if Bethe-Block energy loss used (not e+/e-)
        ONEPASS  = 0x0400,  //Set if the full first pass track fit finished
        TWOPASS  = 0x0800,  //Set if an iteration of the first fit finished
        PRCALSRCH= 0x1000,  //Set if Pat. Rec. used Cal Energy Centroid
        PRBLNSRCH= 0x2000,  //Set if Pat. Rec. used only Track info.
        TOP      = 0x4000,  //Set if track traj. intercepts top tracker plane
        BOTTOM   = 0x8000   //Set if track traj. intercepts first Cal layer

@endverbatim
The definitions should be fairly stable.
<tr><td> Tkr[1/2]TwrEdge  
<td>F<td>   Distance from tower edge of initial point (0 is halfway between the towers, 
             increases towards center of tower) 
<tr><td> Tkr[1/2]PrjTwrEdge  
<td>F<td>   Distance from tower edge of track extrapolated to the layer upstream 
             of the first layer (See Tkr1TwrEdge.) 
<tr><td> Tkr[1/2]DieEdge  
<td>F<td>   Distance from die (wafer) edge of initial point 
             (0 is halfway between the dies, increases toward center of die)  
<tr><td> Tkr[1/2]KalEne  
<td>F<td>   Kalman energy of track 1; this is the energy determined from the multiple scattering 
            along the track (goes like 1/E). Since it is possible to measure a 
            zero scattering angle, which would lead to infinite energy, 
            the minimum measureable angle, which limits the energy to reasonable values  
<tr><td> Tkr[1/2]ConEne  
<td>F<td>   Energy from PatRec energy tool for track 1. 
             The tool computes the total event energy and then partitions it 
             between the first 2 tracks according to their Kalman energies 
             and energy errors  
<tr><td> Tkr[1/2]KalThetaMS  
<td>F<td>   Multiple scattering angle (radians) referenced to first layer. 
             The contributions from all the layers in the track are adjusted 
             for the predicted energy in each layer, and weighted accordingly. 
             So the result is sensitive to the particle type and 
             the chosen energy-loss mechanism. 		
<tr><td> Tkr[1/2][X/Y/Z]Dir  
<td>F<td>   Track [x/y/z] direction cosine  
<tr><td> Tkr[1/2]Phi  
<td>F<td>   Track phi, radians 
            (direction from which particle comes, not particle direction!) 
            range: (0, 2pi)  
<tr><td> Tkr[1/2]Theta  
<td>F<td>   Track theta, radians (direction ditto)  
<tr><td> Tkr[1/2][X/Y/Z]0  
<td>F<td>   Track [x/y/z] position at first hit  
</table>
@subsection best_only Variables that exist only for best track

<table>
<tr><th> Variable <th> Type  <th> Description				                 
<tr><td> TkrRadLength 
<td>F<td>   Radiation lengths traversed by the best track. 
             This is from half-way thru the initial converter to the lowest bi-plane 
             in the tracker, whether or not the track actually gets to the end. 
<tr><td> TkrTwrEdge 
<td>F<td>   The average distance of the best track from the "edge" of each tray, 
             weighted by radiation lengths traversed. 
            (The edge is a plane halfway between the towers. 
<tr><td> TkrTrackLength 
<td>F<td>   Distance between the start of the best track and the grid, along the track axis. 	
<tr><td> Tkr1TwrGap  
<td>F<td>   Length of track in nominal intertower gap, currently set to 18 mm. 
            Can be a small as zero if track exits through bottom of tracker, 
            and as large as the intertower gap, if track crosses to adjacent tower.  
<tr><td> Tkr1ThetaErr  
<td>F<td>   Error on the measurement of theta  
<tr><td> Tkr1PhiErr  
<td>F<td>   Error on the measurement of phi.  
<tr><td> Tkr1ErrAsym  
<td>F<td>   Tkr1SXY/(Tkr1SXX + Tkr1SYY)  
<tr><td> Tkr1CovDet  
<td>F<td>   Determinant of the error matrix, 
             but normalized to remove the dependence on cos(theta)          
<tr><td> Tkr1S[XX/YY]  
<td>F<td>   [x-x/y-y] element of the covariance matrix; square of error on [x/y]  
<tr><td> Tkr1SXY  
<td>F<td>   x-y element of the covariance matrix; covariance  
<tr><td> Tkr1ToTFirst  
<td>F<td>   ToT of first hit on best track 
            (All ToT's are adjusted for pathlength in the measuring and non-measuring 
            directions in the strip, and for the strip width.)  
<tr><td> Tkr1ToTAve  
<td>F<td>   Average ToT for the hits on the best track  
<tr><td> Tkr1ToTTrAve  
<td>F<td>   Average ToT for the hits on the best track, excluding the largest and smallest  
<tr><td> Tkr1ToTAsym  
<td>F<td>   Asymmetry between last two and first two ToT's for the best track  
<tr><td> Tkr1ChisqAsym  
<td>F<td>   Asymmetry between last two and first two track-segment delta-chisquared's  
<tr><td> Tkr1SSDVetoOld  
<td>F<td>   Number of silicon planes between the top of the extrapolated track 
            and the first plane that has a hit near the track. Only planes that have
            wafers which intersect the extrapolated track are considered. No checks 
            for dead strips, etc. are made (yet!).
            Can be used as a back-up for the ACD. 
<tr><td> Tkr1SSDVeto  
<td>F<td>   New version of the SSD Veto. For this variable, tracks which pass close
            to a dead plane, buffer-saturated region, inter-wafer gap, gap between towers,
            or dead strips, do not cause the veto count to be incremented if no hit is found.
            This almost certainly overdoes it: a more correct calculation would include
            the probability for the track to cross the inactive region. (Coming soon!)
<tr><td> TkrVetoPlaneCrossed
<td>I<td>   Number of planes contributiong to the SSD Veto. This doesn't count the points
            where a track crosses in a gap.
<tr><td> Tkr1CoreHC
<td>F<td>   Number of clusters within a roughly cylindrical region )(default radius 10 mm) 
            around the hits in each plane between the first and last on the best
            track, excluding the clusters that belong to the track itself
<tr><td> TkrDispersion
<td>F<td>   the RMS of the distances between the 1st track and all others in the event.
            For tracks which start "above" the first track distance is the 3-D distance
            For the rest, distance is the doca of the head of the track to the axis of
            the first track.
<tr><td> TkrUpstreamHC
<td>F<td>   The number of hits in a cylinder (default radius 150 mm) up to 4 layers thick
            above the head of the first track.
<tr><td> Tkr1CORERatio
<td>F<td>   the ratio of Tkr1CoreHC and Tkr1Hits
<tr><td> Tkr1LATEdge
<td>F<td>   Minimum distance to any LAT edge of the head of the best track
</table>

@subsection ssdveto Diagnostic SSD Veto Variables (Optional, present by default) 
@verbatim
(Turn off with jO parameter: ToolSvc.TkrValsTool.enableVetoDiagnostics = false;)
@endverbatim
     <table>
<tr><th> Variable <th> Type  <th> Description				                 
<tr><td> TkrVetoTrials
<td>I<td>   Difference between the plane number of the last plane crossed and the first, plus one
            Any gaps above the last plane are not counted. (This may change soon.)
<tr><td> TkrVetoHitFound
<td>I<td>   Number of hits found
<tr><td> TkrVetoUnknown
<td>I<td>   Missing hits not ascribable to any inactive area
<tr><td> TkrVetoTower
<td>I<td>   Number of tower crossings (zero for now, because the propagator doesn't report them)
<tr><td> TkrVetoGapCorner
<td>I<td>   Number of missing hits close to a wafer corner
<tr><td> TkrVetoGapEdge
<td>I<td>   Number of missing hits close to a wafer edge
<tr><td> TkrVetoBadCluster
<td>I<td>   Number of missing hits close to a dead strip
<tr><td> TkrVetoDeadPlane
<td>I<td>   Number of missing hits close to a dead plane
<tr><td> TkrVetoTruncated
<td>I<td>   Number of missing hits close to a truncated region
</table>

*/


StatusCode TkrValsTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;

    MsgStream log(msgSvc(), name());

    log << MSG::INFO  << "#################" << endreq << "# ";
    log << (m_useNew ? "New " : "Old ");
    log << "version" << endreq << "#################" << endreq;

    if((ValBase::initialize()).isFailure()) return StatusCode::FAILURE;

    // get the services

    if( serviceLocator()) {

        if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return fail;
        }

        m_towerPitch = m_tkrGeom->towerPitch();
        m_xNum       = m_tkrGeom->numXTowers();
        m_yNum       = m_tkrGeom->numYTowers();
        m_activeWidth = m_tkrGeom->nWaferAcross()*m_tkrGeom->waferPitch() + 
            (m_tkrGeom->nWaferAcross()-1)*m_tkrGeom->ladderGap();

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

    if (toolSvc()->retrieveTool("TkrFlagHitsTool", pFlagHits).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrFlagHitsTool" << endreq;
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
    addItem("TkrSurplusHCInside", &Tkr_SurplusHCInside);
    addItem("TkrSurplusHitRatio", &Tkr_SurplusHitRatio);
    addItem("TkrUpstreamHC",  &Tkr_UpstreamHC);
    addItem("TkrDispersion", &TkrDispersion);

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
    addItem("Tkr1SSDVetoOld", &Tkr_1_SSDVetoOld);
    addItem("Tkr1SSDVeto",    &Tkr_1_SSDVeto);
    addItem("TkrPlaneCrossed",  &Tkr_1_VetoPlaneCrossed);

    if(m_enableVetoDiagnostics) {
        addItem("TkrVetoHitFound",   &Tkr_1_VetoHitFound);
        addItem("TkrVetoTrials",     &Tkr_1_VetoTrials);
        addItem("TkrVetoUnknown",    &Tkr_1_VetoUnknown);
        addItem("TkrVetoTower",      &Tkr_1_VetoTower);
        addItem("TkrVetoGapCorner",  &Tkr_1_VetoGapCorner);
        //addItem("TkrMinGapDistance",  &Tkr_1_MinGapDistance);
        //addItem("TkrMaxGapDistance",  &Tkr_1_MaxGapDistance);
        addItem("TkrVetoGapEdge",    &Tkr_1_VetoGapEdge);
        addItem("TkrVetoBadCluster", &Tkr_1_VetoBadCluster);
        addItem("TkrVetoDeadPlane",  &Tkr_1_VetoDeadPlane);
        addItem("TkrVetoTruncated",  &Tkr_1_VetoTruncated);
    }

    addItem("Tkr1CoreHC",     &Tkr_1_CoreHC);
    addItem("Tkr1LATEdge",    &Tkr_1_LATEdge);

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
    /*
    addItem("TkrFloat", &Tkr_float);
    addItem("TkrInt",   &Tkr_int);
    //try just aliasing it to the first of the three items
    //  is it guaranteed to work?
    addItem("Tkr1Pos[3]",        &Tkr_1_x0);
    addItem("TkrIntArray[5]",    Tkr_ints);
    addItem("TkrDoubleArray[2]", Tkr_doubles);
    addItem("TkrUIntArray[7]",   Tkr_uInts);
    addItem("TkrString",           Tkr_string);
    */

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

    // params for cone model
    double _eConeMin   = 30.;
    double _eConeBreak = 100.;
    double _eConeMax   = 1000.;
    double _coneAngleEMax= 2.65;
    double _coneAngleEBreak = 8.8;
    double _coneAngleEMin = 12.9;
    double _coneOffset    = 2.0;
    double _layerFactor   = 1.0;
    double _expCosth  = 1.5;

    //regions for various hit counts
    double _nearRegion     = 30.0;   // for TkrHDCount
    double _upstreamRegion = 150.0;  // for TkrUpstreamHC
    int    _nUpstream      = 4;      // max number of layers to look upstream
    double _coreRegion     = 10.0;   // for Tkr1CoreHC
    double _vetoRegion     = 10.0;   // for Tkr1SSDVeto
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
        Tkr_1_Hits         = track_1->getNumFitHits();
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

        //test
        /* 
        float fVal;
        getVal("Tkr1Pos[0]", fVal, -1);
        std::cout << "Tkr1Pos[0] " << fVal << std::endl;
        getVal("Tkr1Pos[1]", fVal, -1);
        std::cout << "Tkr1Pos[1] " << fVal << std::endl;
        getVal("Tkr1Pos[2]", fVal, -1);
        std::cout << "Tkr1Pos[2] " << fVal << std::endl;
        //getVal("Tkr1Pos[3]", fVal, -1);                   // illegal argument
        //std::cout << "Tkr1Pos[3] " << fVal << std::endl; 
        //getVal("Tkr1Pos[3", fVal, -1);                    // malformed (missing "]")

        //getVal("TkrFoo", fVal, -1);                       // doesn't exist

        int i;
        for (i=0; i<5;++i) {Tkr_ints[i] = i;}
        for (i=0; i<2;++i) {Tkr_doubles[i] = i+1;}
        for (i=0; i<7;++i) {Tkr_uInts[i] = i+1;}
        

        //std::cout << "xdir = " << Tkr_1_xdir<< ", Before: " << Tkr_string << std::endl;
        std::string tempStr;
        tempStr = ((Tkr_1_xdir>0) ?
            tempStr = "positive num, but way too long for the field" : "neg");
            //Tkr_string[19] = 0;
        strncpy(Tkr_string, tempStr.c_str(), 20);
        Tkr_string[19] = 0;
        //std::cout << "After: " << Tkr_string << std::endl;
        */

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
        Tkr_1_ThetaErr      = t1.z()*t1.z()*sqrt(std::max(0.0, cosPhi*cosPhi*Tkr_1_Sxx + 
            2.*sinPhi*cosPhi*Tkr_1_Sxy + sinPhi*sinPhi*Tkr_1_Syy)); 
        Tkr_1_PhiErr        = (-t1.z())*sqrt(std::max(0.0, sinPhi*sinPhi*Tkr_1_Sxx - 
            2.*sinPhi*cosPhi*Tkr_1_Sxy + cosPhi*cosPhi*Tkr_1_Syy));
        Tkr_1_ErrAsym     = fabs(Tkr_1_Sxy/(Tkr_1_Sxx + Tkr_1_Syy));
        Tkr_1_CovDet      = sqrt(std::max(0.0f,Tkr_1_Sxx*Tkr_1_Syy-Tkr_1_Sxy*Tkr_1_Sxy))*Tkr_1_zdir*Tkr_1_zdir;

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
        double first_ToTs = 0.; 
        double last_ToTs  = 0.;
        double min_ToT   = 1000.; 
        double max_ToT   = -1000.;  
        int    hit_counter = 0; 
        double chisq_first = 0.;
        double chisq_last  = 0.; 

        Event::TkrTrackHitVecConItr pHit = track_1->begin();

        // loop over the hits to calculate various numbers
        double tkrTrackEnergy1 = 0, tkrTrackEnergy2 = 0;
        int plane = m_tkrGeom->getPlane((*pHit)->getTkrId());
        int gapId = -1;
        bool gapFound = false;

        
        // count the number of real hits... may not be necessary but I'm nervous!
        int clustersOnTrack = 0;
        while(pHit != track_1->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            unsigned int bits = hit->getStatusBits();
            if((bits & Event::TkrTrackHit::HITISSSD)==0) continue;
            const Event::TkrCluster* cluster = hit->getClusterPtr();
            double mips = cluster->getMips();
            if (mips<0.0 || mips>10.0) continue;
            clustersOnTrack++;
        }

        pHit = track_1->begin();
        const Event::TkrTrackParams params((*pHit)->getTrackParams(Event::TkrTrackHit::SMOOTHED));
        while(pHit != track_1->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            unsigned int bits = hit->getStatusBits();

            int layer = m_tkrGeom->getLayer(hit->getTkrId());
            convType type = m_tkrGeom->getLayerType(layer);
            if (type==STANDARD)   {tkrTrackEnergy1 += cfThin;}
            else if (type==SUPER) {tkrTrackEnergy1 += cfThick;}

            // check if hit is in an ssd
            if ( !gapFound && (bits & Event::TkrTrackHit::HITISSSD)==0) {
                Point  gapPos = hit->getPoint(Event::TkrTrackHit::PREDICTED);
                Tkr_1_GapX = gapPos.x();
                Tkr_1_GapY = gapPos.y();
                idents::TkrId tkrId = hit->getTkrId();
                if(tkrId.hasTray()) {
                    gapId = m_tkrGeom->getPlane(hit->getTkrId());
                } else {
                    gapId = plane;
                }
                gapFound = true;
            } else {
            }
            plane--;
            if ((bits & Event::TkrTrackHit::HITISSSD)==0) continue;
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
            if(mips<0.0 || mips>10.0) continue;

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
            // first 2 valid clusters
            if(hit_counter < 3) {
                first_ToTs += mips;
                chisq_first += hit->getChiSquareSmooth();
            }
            // last 2 valid clusters
            if(hit_counter > clustersOnTrack-2){
                last_ToTs += mips;
                chisq_last += hit->getChiSquareSmooth();
            }
        }

        if(clustersOnTrack>3) {
            Tkr_1_ToTTrAve = (Tkr_1_ToTAve - max_ToT - min_ToT)/(clustersOnTrack-2.);
            Tkr_1_ToTAve /= clustersOnTrack;
            if(first_ToTs+last_ToTs>0) {
                Tkr_1_ToTAsym = (last_ToTs - first_ToTs)/(first_ToTs + last_ToTs);
            }
        }

        tkrTrackEnergy1 /= fabs(Tkr_1_zdir);


        Tkr_1_FirstGapPlane = gapId; 

        // Chisq Asymmetry - Front vs Back ends of tracks
        if (chisq_last+chisq_first>0) Tkr_1_ChisqAsym = (chisq_last - chisq_first)/(chisq_last + chisq_first);

        //m_G4PropTool->setStepStart(x1, -t1); //Note minus sign - swim backwards towards ACD
  
        bool upwards = true;
        m_G4PropTool->setStepStart(params, x1.z(), upwards);

        int topPlane = m_tkrGeom->numPlanes()-1; 
        double topOfTkr = m_tkrGeom->getPlaneZ(topPlane) + 1.0;
        double arc_min = fabs((topOfTkr-x1.z())/t1.z());
        arc_min = std::min( arc_min, maxPath); 
        m_G4PropTool->step(arc_min);  
        int numSteps = m_G4PropTool->getNumberSteps();

        idents::VolumeIdentifier volId;
        idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();

        int firstPlane = m_tkrGeom->getPlane(track_1->front()->getTkrId()); 
        int firstLayer = m_tkrGeom->getLayer(firstPlane);
        double zFirstLayer = m_tkrGeom->getLayerZ(firstLayer);

        //int numLayersTrack = (int) (Tkr_1_FirstLayer - Tkr_1_LastLayer) + 1;

        double costh = fabs(t1.z());
        double secth = 1./costh;

        // for the footprints
        double secthX = 1./sqrt(1.0 - t1.x()*t1.x());
        double secthY = 1./sqrt(1.0 - t1.y()*t1.y());

        double xVetoRgn = _vetoRegion*secthX;
        double yVetoRgn = _vetoRegion*secthY;

        // Note: skip the first vol - as this is the head SSD

        double arcLen = m_G4PropTool->getStepArcLen(0);
        bool isFirstPlane = true;

     try {
        for(int istep = 1; istep < numSteps; ++istep) { 
            volId = m_G4PropTool->getStepVolumeId(istep);
            volId.prepend(prefix);
            Point x_step       = m_G4PropTool->getStepPosition(istep);
            arcLen += m_G4PropTool->getStepArcLen(istep);

            bool forward = false;
            const Event::TkrTrackParams newParams = 
                m_G4PropTool->getTrackParams(arcLen, Tkr_1_ConEne, forward);
 
            //std::cout << "pos " << x_step << ", step " << istep << ", arcLen " << arcLen << std::endl;
            //std::cout << "err " << sqrt(newParams(xPosIdx,xPosIdx)) << " "
            //        << sqrt(newParams(yPosIdx,yPosIdx)) << std::endl;
            // we're outside the LAT
            if(x_step.z() > topOfTkr || !m_tkrGeom->isInActiveLAT(x_step) ) break; 

            // check that it's really a TKR hit (probably overkill)
            if(volId.size() != 9) continue; 
            if(!(volId[0]==0 && volId[3]==1)) continue; // !(LAT && TKR)
            if(volId[6]> 1) continue;  //It's a converter or some other tray element!
            Tkr_1_VetoPlaneCrossed++;

            // now check if there's a hit near the extrapolated track!
            // if there is, reset the veto counter... we want leading non-hits
            //std::cout << "Tkr1SSDVeto volId " << volId.name() << std::endl;

            int tower = idents::TowerId(volId[2], volId[1]).id();
            int tray = volId[4];
            int view = volId[5];
            int face = volId[6];
            int layer = m_tkrGeom->trayToBiLayer(tray, face);
            int plane = m_tkrGeom->trayToPlane(tray, face);

            int firstPlane;
            if(isFirstPlane) {
                isFirstPlane = false;
                firstPlane = plane;
            }

            // I think we want to do this, most likely this is missed for
            //   some good reason.
            // on the other hand, it is a missed hit, so remove test for the new code
            
            
            Tkr_1_VetoTrials = abs(plane-firstPlane) + 1;

            //if(!m_useNew&&layer==firstLayer) continue;

            double vetoRgn = (view==0 ? xVetoRgn : yVetoRgn);

            int nVetoHits = pQueryClusters->numberOfHitsNear(view, layer, 
                vetoRgn, x_step, t1);

            if (nVetoHits>0) {
                // found a hit, reset the SSDVeto
                    Tkr_1_SSDVeto = 0.0;
                    Tkr_1_SSDVetoOld = 0.0;
                    if(m_enableVetoDiagnostics) { Tkr_1_VetoHitFound++; }

            } else { 
                // no hit
                // plain old ssdveto
                Tkr_1_SSDVetoOld += 1.0;

                idents::TkrId tkrId(volId);
                Event::TkrTrackParams outParams;

                unsigned int status_bits = 0;
                int stage = -1;
                if(m_useNew||m_enableVetoDiagnostics) {
                    stage = pFlagHits->flagHits(tkrId, newParams, x_step.z(),
                        m_minVetoError, m_maxVetoError, m_vetoNSigma, 
                        outParams, status_bits);
                }
                if(m_useNew && (stage==-1)) { 
                    Tkr_1_SSDVeto += 1.0; 
                }

                if(m_enableVetoDiagnostics) {
                    switch (stage) {
                        case -1:
                            Tkr_1_VetoUnknown++;
                            break;
                        case 1: 
                            Tkr_1_VetoDeadPlane++;
                            break;
                        case 2:
                            Tkr_1_VetoTruncated++;
                            break;
                        case 3:
                            Tkr_1_VetoTower;
                            break;
                        case 4:
                            Tkr_1_VetoGapCorner++;
                            break;
                        case 5:
                            Tkr_1_VetoGapEdge++;
                            break;
                        case 6:
                            Tkr_1_VetoBadCluster++;
                            break;
                        default:
                            std::cout << "shouldn't get here, stage = " 
                                << stage << std::endl;
                    } // end switch
                } // end diagnostics
            } // end no hit 
        } // end loop over steps

    } catch( std::exception& e) {
      SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);
      unsigned long evtId = (header) ? header->event() : 0;
      long runId = (header) ? header->run() : -1;
      log << MSG::WARNING << "Caught exception (run,event): ( "
          << runId << ", " << evtId << " ) " << e.what()
          << " Skipping the TKR calculations" << endreq;

      Tkr_1_VetoPlaneCrossed = s_badVal; 
      Tkr_1_VetoTrials = s_badVal;
      Tkr_1_SSDVeto = s_badVal;
      Tkr_1_SSDVetoOld = s_badVal;
      Tkr_1_VetoUnknown = s_badVal;
      Tkr_1_VetoDeadPlane = s_badVal;
      Tkr_1_VetoTruncated = s_badVal; 
      Tkr_1_VetoTower = s_badVal;   
      Tkr_1_VetoGapCorner = s_badVal;
      Tkr_1_VetoGapEdge = s_badVal;   
      Tkr_1_VetoBadCluster = s_badVal;

    } catch (...) {
      SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);
      unsigned long evtId = (header) ? header->event() : 0;
      long runId = (header) ? header->run() : -1;
      log << MSG::WARNING << "Caught exception (run,event): ( "
          << runId << ", " << evtId << " ) " 
          << " Skipping the TKR calculations" << endreq;
      Tkr_1_VetoPlaneCrossed = s_badVal; 
      Tkr_1_VetoTrials = s_badVal;
      Tkr_1_SSDVeto = s_badVal;
      Tkr_1_SSDVetoOld = s_badVal;
      Tkr_1_VetoUnknown = s_badVal;
      Tkr_1_VetoDeadPlane = s_badVal;
      Tkr_1_VetoTruncated = s_badVal; 
      Tkr_1_VetoTower = s_badVal;   
      Tkr_1_VetoGapCorner = s_badVal;
      Tkr_1_VetoGapEdge = s_badVal;   
      Tkr_1_VetoBadCluster = s_badVal;
    }

        // minimum distance from any edge, measured from the edge of the active area
        double deltaEdge = 0.5*(m_towerPitch - m_tkrGeom->trayWidth()) 
            - m_tkrGeom->siDeadDistance() ;
        double tkrXLo = m_tkrGeom->getLATLimit(0,LOW)  + deltaEdge;
        double tkrXHi = m_tkrGeom->getLATLimit(0,HIGH) - deltaEdge;
        double tkrYLo = m_tkrGeom->getLATLimit(1,LOW)  + deltaEdge;
        double tkrYHi = m_tkrGeom->getLATLimit(1,HIGH) - deltaEdge;

        double xEdge = std::min(Tkr_1_x0-tkrXLo, tkrXHi- Tkr_1_x0);
        double yEdge = std::min(Tkr_1_y0-tkrYLo, tkrYHi- Tkr_1_y0);

        Tkr_1_LATEdge = (float) std::min(xEdge, yEdge);

        if(nTracks > 1) {
            pTrack++;

            // try Bill's dispersion variable here

            Doca trk1Doca(x1, t1);

            Event::TkrTrackColConPtr trkIter;

            for(trkIter=pTrack; trkIter!=pTracks->end(); ++trkIter) {
                Event::TkrTrack* trk = *trkIter;
                double docaTrk = trk1Doca.docaOfPoint(trk->getInitialPosition());
                TkrDispersion += (float) docaTrk*docaTrk;
                double s = trk1Doca.arcLenRay1();
                if (s<0) {
                    TkrDispersion += (float) s*s;
                }
            }
            TkrDispersion = sqrt(TkrDispersion/(nTracks-1));


            const Event::TkrTrack* track_2 = *pTrack;
            Tkr_2_Chisq        = track_2->getChiSquareSmooth();
            Tkr_2_FirstChisq     = track_2->chiSquareSegment();
            Tkr_2_FirstGaps      = track_2->getNumXFirstGaps() + track_2->getNumYFirstGaps();
            Tkr_2_Qual         = track_2->getQuality();
            Tkr_2_Type         = track_2->getStatusBits();
            Tkr_2_Hits         = track_2->getNumFitHits();
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

            Event::TkrTrackHitVecConItr pHit = track_2->begin();

            // count up the hits for track energy
            while(pHit != track_2->end()) {
                const Event::TkrTrackHit* hit = *pHit++;

                int layer = m_tkrGeom->getLayer(hit->getTkrId());
                convType type = m_tkrGeom->getLayerType(layer);
                if (type==STANDARD)   {tkrTrackEnergy2 += cfThin;}
                else if (type==SUPER) {tkrTrackEnergy2 += cfThick;}
            }
            tkrTrackEnergy2 /= fabs(Tkr_2_zdir);
        }

        Tkr_Sum_KalEne    = Tkr_1_KalEne+Tkr_2_KalEne; 
        Tkr_Sum_ConEne    = Tkr_1_ConEne+Tkr_2_ConEne;      

        double tkrTrackEnergy = tkrTrackEnergy1 + tkrTrackEnergy2;

        // Computation of the tracker contribution to the total energy 
        arc_min = (x1.z() - m_tkrGeom->calZTop())*secth; 
        m_G4PropTool->setStepStart(x1, t1);
        m_G4PropTool->step(arc_min);
        double z_present = x1.z();

        // Compute the sum-of radiation_lengths x Hits in each layer
        //double tracker_ene_corr = 0.; 
        double rad_len_sum  = 0.; 
        double radlen       = 0.;
        double radlen_old   = 0.; 
        double arc_len      = 0.; 
        int    total_hits   = 0; 
        int    thin_hits    = 0;
        int    thick_hits   = 0; 
        int    blank_hits   = 0; 
        double ave_edge     = 0.; 

        //float surplus_in = 0;
        //float total_layer_hits = 0;
        int numTowers = m_xNum*m_yNum;
        std::vector<float> layerInCount(numTowers,0.0);
        std::vector<float> layerOutCount(numTowers,0.0);

        // do the hit counts
        // Tkr_HDCount

        double xNearRgn = _nearRegion*secthX;
        double yNearRgn = _nearRegion*secthY;

        Tkr_HDCount = pQueryClusters->numberOfUUHitsNear((int) Tkr_1_FirstLayer, 
            xNearRgn, yNearRgn, x1);
        
        // Tkr1CoreHC:

        double xCoreRgn = _coreRegion*secthX;
        double yCoreRgn = _coreRegion*secthY;

        int numLayers = m_tkrGeom->numLayers();

        Event::TkrTrackHitVecConItr hitIter = track_1->begin();
        for(;hitIter!=track_1->end(); ++hitIter) {
            const Event::TkrTrackHit* thisHit = *hitIter;
            idents::TkrId thisId = thisHit->getTkrId();
            int thisPlane = m_tkrGeom->getPlane(thisId);
            int thisLayer = m_tkrGeom->getLayer(thisPlane);
            int thisView  = thisId.getView();
            Point pos;
            if(thisHit->validMeasuredHit()) {
                pos  = thisHit->getPoint(Event::TkrTrackHit::MEASURED);
            } else if(thisHit->validFilteredHit()) {
                pos  = thisHit->getPoint(Event::TkrTrackHit::FILTERED);
            } else { continue; }

            double distance = (thisView==0 ? xCoreRgn : yCoreRgn);
            int coreHits = pQueryClusters->numberOfHitsNear(thisView, thisLayer,
                distance, pos);
            if (thisHit->validCluster()) coreHits--;
            Tkr_1_CoreHC += coreHits;
        }

        //TkrUpstreamHC

        int topLayer  = numLayers - 1;
        int layer    = firstLayer+1;
        int upperLayer = std::min(firstLayer+_nUpstream, topLayer);
        double xUpstreamRgn = _upstreamRegion*secthX;
        double yUpstreamRgn = _upstreamRegion*secthY;

        for(; layer<=upperLayer; ++layer) {
            double zLayer = m_tkrGeom->getLayerZ(layer);
            double deltaZ = zFirstLayer - zLayer;
            double arcLength = deltaZ*secth;
            Point x_hit = x1 + arcLength*t1;
            Tkr_UpstreamHC += pQueryClusters->numberOfHitsNear(layer, 
                xUpstreamRgn, yUpstreamRgn, x_hit, t1);       
        }

        // Surplus hits

        double tanth = (1.0-costh*costh)*secth;
        double spread0 = _coneOffset + _layerFactor*tanth;
        double cosFactor = pow(secth, _expCosth);

        float Tkr_SurplusHCOutside = 0.0f;
        //float Tkr_SurplusHCInside  = 0.0f;
        //float Tkr_Total_Hits = 0.0f;

        // hate to do this, but we need ERecon
        // Recover pointer to CalEventEnergy info 
        double CAL_EnergyCorr = 0.0;
#ifdef PRE_CALMOD
       Event::CalEventEnergy* calEventEnergy = 
            SmartDataPtr<Event::CalEventEnergy>(m_pEventSvc, EventModel::CalRecon::CalEventEnergy);
#else
        Event::CalEventEnergyCol * calEventEnergyCol = 
          SmartDataPtr<Event::CalEventEnergyCol>(m_pEventSvc, EventModel::CalRecon::CalEventEnergyCol);
        Event::CalEventEnergy * calEventEnergy = 0 ;
        if ((calEventEnergyCol!=0)&&(!calEventEnergyCol->empty()))
          calEventEnergy = calEventEnergyCol->front() ;
#endif
        if (calEventEnergy != 0) {
            // Extraction of results from CalValCorrTool in CalRecon... 
            Event::CalCorToolResultCol::iterator corIter = calEventEnergy->begin();
            for(;corIter != calEventEnergy->end(); corIter++){
                Event::CalCorToolResult corResult = **corIter;
                if (corResult.getCorrectionName() == "CalValsCorrTool") {
                    CAL_EnergyCorr   = corResult["CorrectedEnergy"];
                }
            }
        }
        
        double eRecon = CAL_EnergyCorr + tkrTrackEnergy;
        double eCone = std::min(_eConeMax, std::max(eRecon, _eConeMin));
        // Get the basic cone angle for this Energy
        double coneAngle;
        if (eCone<_eConeBreak) {
            coneAngle = interpolate(1.0/eCone, 1.0/_eConeMin, 1.0/_eConeBreak, 
            _coneAngleEMin, _coneAngleEBreak);
        } else {
            coneAngle = interpolate(1./eCone, 1./_eConeBreak, 1.0/_eConeMax, 
            _coneAngleEBreak, _coneAngleEMax);
        }

        double xSprd0 = coneAngle*secthX*cosFactor;
        double ySprd0 = coneAngle*secthY*cosFactor;

        for(layer = firstLayer; layer>=0; --layer) {
            
            if(layer <firstLayer) {
                radlen = m_G4PropTool->getRadLength(arc_len); 
            }

            // Assume location of shower center is given by 1st track

            // try to get actual x and y (accounting for the differences in z)
            double zX = m_tkrGeom->getLayerZ(layer, 0);
            double zY = m_tkrGeom->getLayerZ(layer, 1);
            double zAve = 0.5*(zX + zY);
            double arcLenX = (x1.z() - zX)*secth;
            double arcLenY = (x1.z() - zY)*secth;

            double x_hitX = x1.x() + arcLenX*t1.x();
            double x_hitY = x1.y() + arcLenY*t1.y();
            Point x_hit1(x_hitX, x_hitY, zAve);
            // whew!!

            int hitXTower, hitYTower;
            m_tkrGeom->truncateCoord(x_hit1.x(), m_towerPitch, m_xNum, hitXTower);
            m_tkrGeom->truncateCoord(x_hit1.y(), m_towerPitch, m_yNum, hitYTower);
            int thisTower = idents::TowerId(hitXTower, hitYTower).id();
                
            Point x_hit = x1 + arc_len*t1;

            // trial code for Surplus hits

            layerInCount.assign(numTowers,0);
            layerOutCount.assign(numTowers, 0);

            Event::TkrClusterVec clusVec[2];
            int tower;
            std::vector<int> clusCount[2];
            clusCount[0].assign(numTowers,0);
            clusCount[1].assign(numTowers,0);
            int view;
            // fpr each layer, count the clusters in each tower, view
            for (view=0; view<2; ++view) {
                clusVec[view] = pQueryClusters->getClusters(view, layer);
                //std::cout << clusVec[view].size() << std::endl;
                Event::TkrClusterVecConItr iter = clusVec[view].begin();
                for(;iter!=clusVec[view].end(); ++iter) {
                    idents::TkrId id = (*iter)->getTkrId();
                    tower = idents::TowerId(id.getTowerX(), id.getTowerY()).id();
                    clusCount[view][tower]++;
                    //std::cout << "count, " << tower << " " << view << ": " 
                    //    << clusCount[view][tower] << std::endl;
                }
            }

            // form the x-y coincidence
            std::vector<bool> isXY(numTowers, false);
            for (tower=0;tower<numTowers; ++ tower) {
                isXY[tower] = (clusCount[0][tower]>0 && clusCount[1][tower]>0);
                //std::cout << "tower " << tower << ", " << clusCount[0][tower] 
                //    << " " << clusCount[1][tower] << ", " << isXY[tower] << std::endl;
            }
            // make sure *this* tower is included!
            isXY[thisTower] = true;

            float factor;

            double xSprd = spread0 + xSprd0*(firstLayer-layer);
            double ySprd = spread0 + ySprd0*(firstLayer-layer);

            // test each cluster in surviving towers
            int mode = 0;
            if (mode==0) {
                factor = 1.0;
                for (view=0; view<2; ++view) {
                    Event::TkrClusterVecConItr iter = clusVec[view].begin();
                    for(;iter!=clusVec[view].end(); ++iter) {
                        idents::TkrId id = (*iter)->getTkrId();
                        tower = idents::TowerId(id.getTowerX(), id.getTowerY()).id();
                        if(!isXY[tower]) continue;
                        Vector diff = x_hit1 - (*iter)->position();
                        bool in;
                        // replace with current definition of "in"
                        if(view==0) {
                            in = (fabs(diff.x())<=xSprd && (fabs(diff.y())-ySprd<0.5*m_activeWidth));
                        } else { 
                            in = (fabs(diff.y())<=ySprd && (fabs(diff.x())-xSprd<0.5*m_activeWidth));
                        }
                        if (in) {layerInCount[tower]  += 1.0f;}
                        else    {layerOutCount[tower] += 1.0f;}
                    }
                }
            } 
            /*
            // first attempt at code that uses TkrPoint-like objects... needs lots more work!!
            else if (mode==1) {
                double xDenom = 1./xSprd/xSprd;
                double yDenom = 1./ySprd/ySprd;
                Event::TkrClusterVecConItr iterX = clusVec[0].begin();
                for(;iterX!=clusVec[0].end(); ++iterX) {
                    idents::TkrId idX = (*iterX)->getTkrId();
                    int towerX = idents::TowerId(idX.getTowerX(), idX.getTowerY()).id();
                    if(!isXY[tower]) continue;
                    Event::TkrClusterVecConItr iterY = clusVec[1].begin();
                    for(;iterY!=clusVec[1].end(); ++iterY) {
                        idents::TkrId idY = (*iterY)->getTkrId();
                        int towerY = idents::TowerId(idY.getTowerX(), idY.getTowerY()).id();
                        if(towerX!=towerY) continue;
                        layerOutCount[tower] += 1.0f;
                        double dx = fabs(x_hit.x() - (*iterX)->position().x());
                        double dy = fabs(x_hit.y() - (*iterY)->position().y());
                        if (dx>xSprd || dy>ySprd) continue;
                        // could be inside!
                        if(dx*dx/xDenom + dy*dy/yDenom < 1.) layerInCount[tower] += 1.0f;
                    }
                }
                //    factor = ((float)clusCount[0][tower]+clusCount[1][tower])/
                //    std::max(clusCount[0][tower]*clusCount[1][tower], 1);
            } */        

            int numHits = 0, numHitsOut = 0;
            for(tower=0;tower<numTowers;++tower) {
                numHitsOut += (int)layerOutCount[tower];
                numHits    += (int)(layerInCount[tower]*factor);
            }
            Tkr_SurplusHCOutside += numHitsOut;
            Tkr_SurplusHCInside  += numHits;

            double layer_edge = towerEdge(x_hit);

            double delta_rad= radlen-radlen_old;
            double thisRad;
            //A bit cleaner
            switch (m_tkrGeom->getLayerType(layer)) {
                case STANDARD:
                    thisRad = radThin;
                    thin_hits += numHits;
                    break;
                case SUPER:
                    thisRad = radThick;
                    thick_hits += numHits;
                    break;
                case NOCONV:
                    thisRad = 0.0;
                    blank_hits += numHits;
                    break;
                default:
                    break;
            }
            if (layer==firstLayer) {
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
            if(layer==0) break;

            double z_next = m_tkrGeom->getLayerZ(layer-1);
            double deltaZ = z_present - z_next;
            z_present = z_next;

            arc_len += fabs( deltaZ/t1.z()); 
            radlen_old = radlen; 
        }

        Tkr_Energy     = (cfThin*thin_hits + cfThick*thick_hits)*std::min(secth, 5.0);
        Tkr_SurplusHitRatio = Tkr_SurplusHCOutside/std::max(1.0f, Tkr_SurplusHCInside);        

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
