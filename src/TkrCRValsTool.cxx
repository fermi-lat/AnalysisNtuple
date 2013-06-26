/** @file TkrCRValsTool.cxx
@brief Calculates the Tkr analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrCRValsTool.cxx,v 1.2 2013/05/16 23:44:16 lsrea Exp $
*/
//#define PRE_CALMOD 1

// To Do:
// implement better code to check if in tower

// Include files



#include "AnalysisNtuple/IValsTool.h"
#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataSvc.h"
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

/*! @class TkrCRValsTool
@brief calculates Cosmic Ray Finder Tracker values

@authors Bill Atwood, Leon Rochester
*/

class TkrCRValsTool :  public ValBase 
{
public:

    TkrCRValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~TkrCRValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

private:

    double towerEdge(Point pos) const;
    double containedFraction(Point pos, double gap, double r, 
        double costh, double phi) const;
    float SSDEvaluation(const Event::TkrTrack* track); 

    // some local constants
    double m_towerPitch;
    int    m_xNum;
    int    m_yNum;
    double m_activeWidth;
    bool   m_useNew;
    bool   m_enableVetoDiagnostics;
    int    m_messageCount;

    // Local variables to transfer results of SSD calculation
    float m_VetoPlaneCrossed; 
    float m_VetoTrials;
    float m_SSDVeto;
    float m_VetoUnknown;
    float m_VetoDeadPlane;
    float m_VetoTruncated; 
    float m_VetoTower;   
    float m_VetoGapCorner;
    float m_VetoGapEdge;   
    float m_VetoBadCluster;
    float m_VetoHitFound;

    // some pointers to services
	/// pointer to data service
	DataSvc*         m_dataSvc;

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

    IPropagator* m_G4PropTool; 
    /// AcdValsTool for Veto Track Number
    IValsTool* m_pAcdTool;

    // properties
    double m_minVetoError;
    double m_maxVetoError;
    double m_vetoNSigma;
    bool   m_testExceptions;
    int    m_minWide;
    int    m_minWider;


    //Global Track Tuple Items
    float TkrCR_Num_Tracks;
    float TkrCR_TrackLength;

    //First Track Specifics
    float TkrCR_1_Chisq;
    float TkrCR_1_FirstChisq;
    float TkrCR_1_Gaps;
    float TkrCR_1_FirstGapPlane; 
    float TkrCR_1_GapX;
    float TkrCR_1_GapY;
    float TkrCR_1_FirstGaps; 
    float TkrCR_1_Hits;
    float TkrCR_1_FirstHits;
    float TkrCR_1_FirstLayer; 
    float TkrCR_1_LastLayer;

    float TkrCR_1_SaturatedFrac;
    float TkrCR_1_ToT255Frac;
    float TkrCR_1_BothFrac;
    float TkrCR_1_GhostFrac;
    float TkrCR_1_DiagFrac;
    float TkrCR_1_WideFrac;
    float TkrCR_1_WiderFrac;

    float TkrCR_1_Qual;
    float TkrCR_1_Type;

    float TkrCR_1_DifHits;
    float TkrCR_1_KalEne;
    float TkrCR_1_ConEne;
    float TkrCR_1_KalThetaMS;
    float TkrCR_1_RangeEne;
    float TkrCR_1_TwrEdge;
    float TkrCR_1_PrjTwrEdge;
    float TkrCR_1_DieEdge;
    float TkrCR_1_TwrGap;
    double TkrCR_1_xdir;
    double TkrCR_1_ydir;
    double TkrCR_1_zdir;
    float TkrCR_1_Phi;
    float TkrCR_1_Theta;
    float TkrCR_1_x0;
    float TkrCR_1_y0;
    float TkrCR_1_z0;
 
    float TkrCR_1_ToTFirst;
    float TkrCR_1_ToTAve;
    float TkrCR_1_ToTTrAve;
    float TkrCR_1_ToTAsym;
    float TkrCR_1_ChisqAsym;
    float TkrCR_1_SSDVeto; 

    // for SSDVeto Diagnostics
    int   TkrCR_1_VetoTrials;
    int   TkrCR_1_VetoHitFound;
    int   TkrCR_1_VetoUnknown;
    int   TkrCR_1_VetoPlaneCrossed;
    int   TkrCR_1_VetoTower;
    int   TkrCR_1_VetoGapCorner;
    //double Tkr_1_MinGapDistance;
    //double Tkr_1_MaxGapDistance;
    int   TkrCR_1_VetoGapEdge;
    int   TkrCR_1_VetoBadCluster;
    int   TkrCR_1_VetoDeadPlane;
    int   TkrCR_1_VetoTruncated;

    float TkrCR_1_CoreHC;
    float TkrCR_1_LATEdge;

    float TkrCR_Veto_SSDVeto;

};

namespace 
{
    double interpolate ( double xi, double xmin, double xmax, 
        double ymin, double ymax) 
    {
        double step = (xi-xmin)/(xmax-xmin);
        return ymin + (ymax-ymin)*step;
    }

    int xPosIdx = Event::TkrTrackParams::xPosIdx;
    int yPosIdx = Event::TkrTrackParams::yPosIdx;
    int xSlpIdx = Event::TkrTrackParams::xSlpIdx;
    int ySlpIdx = Event::TkrTrackParams::ySlpIdx;

    const int    _badInt  =   -1;
    const float _badFloat = -2.0; 
}

// Static factory for instantiation of algtool objects
//static ToolFactory<TkrValsTool> s_factory;
//const IToolFactory& TkrValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrCRValsTool);

// Standard Constructor
TkrCRValsTool::TkrCRValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this);

    //declareProperty("useNew", m_useNew=true);
    // Old SSDVeto variable has been removed; always use "new" code
    m_useNew = true;
    declareProperty("enableVetoDiagnostics", m_enableVetoDiagnostics=false);
    declareProperty("minVetoError", m_minVetoError=1.0);
    declareProperty("maxVetoError", m_maxVetoError=100000.0);
    declareProperty("vetoNSigma", m_vetoNSigma=2.0);
    declareProperty("testExceptions", m_testExceptions=false);
    declareProperty("minWide", m_minWide=5);
    declareProperty("minWider", m_minWider=9);

}

/** @page anatup_vars 
@section tkrvalstool TkrValsTool Variables

Notes: 
- Variables called Tkr1Xxx refer to the "best" track; 
those called Tkr2Xxx refer to the second track.
- A number of variables have the word "Hits" in their name. This <em>usually</em>
refers to clusters! In some cases it refers to TkrTrackHits. 
The description should make it clear which meaning is intended.
- For variables listed as Tkr[1/2]Xxx there are two versions in the ntuple, 
one for the best and one for the second track. 
- The labels are not entirely consistent, but it's probably 
too disruptive to fix them at this point.
For example: TkrRadLength, TkrTrackLength, TkrTwrEdge refer to track 1. 
Also, Tkr2Angle and Tkr2HDoca are quantities that depend on both tracks.
- The variables associated with the second track are undefined 
if there is only one track! 
Check TkrNumTracks before using these variables! 
In fact check TkrNumTracks before using first-track variables, 
for the same reason.
- A new section of (optional) ssd-veto diagnostic variables has been added. 
They are not written out by default.
- Several new variables starting with "TkrV" have been added. These refer to 
quantities associated with the track likely to have cause the Acd veto.
- Some deleted variables, all Tkr2: FirstHits, DifHits, Gaps, FirstGaps,
DieEdge, KalThetaMs, [X/Y/Z]Dir, Phi, Theta, [X/Y/Z]0.

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
<td>F<td>   Tracker edge correction. This may go away; 
it's an intermediate quantity 
<tr><td> TkrHDCount 
<td>F<td>   Number of unused clusters in top x-y layer of the best track 
within a radius of 30 mm, corrected for track angle
(Used in PSF analysis and background rejection) 
<tr><td> TkrTotalHits 
<td>F<td>   Deprecated. Use TkrSurplusHCInside instead
<tr><td> TkrSurplusHitsInside 
<td>F<td>   Number of clusters inside an energy- and angle-dependent cone 
centered on the reconstructed axis of the best track and
starting at the head of track 1. Only clusters in layers with at
least one x and one y cluster in the tower are counted.
<tr><td> TkrSurplusHitRatio
<td>F<td>   Ratio of the number of clusters outside the cone to the number
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
<tr><td> Tkr[1/2]ToT255Frac
<td>F<td>   Fraction of good hits on track that have ToT==255
<tr><td> Tkr[1/2]BothFrac
<td>F<td>   Fraction of good hits on track that are both regular and diagnostic ghosts
<tr><td> Tkr[1/2]GhostFrac
<td>F<td>   Fraction of good hits on track that are regular ghosts
<tr><td> Tkr[1/2]DiagFrac
<td>F<td>   Fraction of good hits on track that are diagnostic ghosts
<tr><td> Tkr[1/2]SaturatedFrac
<td>F<td>   Fraction of good hits on track that are saturated
<tr><td> Tkr[1/2]WideFrac
<td>F<td>   Fraction of good hits on track whose with is >= m_minWide (default=5)
<tr><td> Tkr[1/2]WiderFrac
<td>F<td>   Fraction of good hits on track whose with is >= m_minWider (default=9)
<tr><td> Tkr[1/2][First/Last]Layer  
<td>F<td>   [First/Last] layer in track  (layer 0 is the bottom of the tracker)
<tr><td> Tkr1FirstGapPlane  
<td>F<td>   plane number of first gap on track 1  
(This and the following X,Y pair can be used to find dead strips)
<tr><td> Tkr1[X/Y]Gap  
<td>F<td>   [x/y] location of first gap on track 1  
<tr><td> Tkr[1/2]Qual  
<td>F<td>   Track "quality": depends on the number of clusters and chisquared of the track. 
Maximum is currently 64, can be negative if chisqared gets large. 
This is used primarily to order the tracks during patrec. 
<strong>It's not a good idea to cut on this variable!</strong>  
<tr><td> Tkr[1/2]Type  
<td>F<td>   These are the status bits from the trackign, containing information
about how the track was found and fitted.   
See Event/Recon/TkrRecon/TkrTrack.h for the current description.     
As of GlastRelease v17r31 (June 2009), the status word bits organized as follows:

@verbatim

Low-order bits:
 |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
  [ Pat Rec Info  ] [Pass ] [ E-Loss] [ Track Energy ]  [Track Fit Status]
High-order bits:
 |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
                                                                 [Ghosts]

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

GHOST      = 0x10000, // set if track contains ghost clusters
DIAGNOSTIC = 0x20000  // set if track contains diagnostic ghost clusters
 
@endverbatim
The definitions should be fairly stable.
<tr><td> Tkr[1/2]TwrEdge  
<td>F<td>   Distance from tower edge of initial point (0 is halfway between the towers, 
increases towards center of tower) 
<tr><td> Tkr[1/2]PrjTwrEdge  
<td>F<td>   Distance from tower edge of track extrapolated to the layer upstream 
of the first layer (See Tkr1TwrEdge.) 
<tr><td> Tkr[1/2]KalEne  
<td>F<td>   Kalman energy of track 1; this is the energy determined from the multiple scattering 
along the track (goes like 1/E). Since it is possible to measure a 
zero scattering angle, which would lead to infinite energy, 
the minimum measureable angle, which limits the energy to reasonable values
<tr><td> Tkr[1/2]RangeEnergy
<td>F<td>  Energy estimated for "stopping tracks", which means (so far) a track that ends before
exiting the bottom, with a wide cluster in the last layer
<tr><td> Tkr[1/2]ConEne  
<td>F<td>   Energy from PatRec energy tool for track 1. 
The tool computes the total event energy and then partitions it 
between the first 2 tracks according to their Kalman energies 
and energy errors  
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
<td>F<td>   ToT associated with the first hit on best track 
(All ToT's are adjusted for pathlength in the measuring and non-measuring 
directions in the strip, and for the strip width.) <em>Note: there is only 
one ToT per half-plane. If there is more than one hit strip, the highest
ToT is stored.</em>
<tr><td> Tkr1ToTAve  
<td>F<td>   Average ToT for the hits on the best track (See note above.) 
<tr><td> Tkr1ToTTrAve  
<td>F<td>   Average ToT for the hits on the best track, 
excluding the largest and smallest (See note above.)
<tr><td> Tkr1ToTAsym  
<td>F<td>   Asymmetry between last two and first two ToT's for the best track
(See note above.)
<tr><td> Tkr1ChisqAsym  
<td>F<td>   Asymmetry between last two and first two track-segment delta-chisquared's  
<tr><td> Tkr1SSDVetoOld 
<td>F<td>   Number of silicon planes between the top of the extrapolated track 
and the first plane that has a cluster near the track. Only planes that have
wafers which intersect the extrapolated track are considered. No checks 
for dead strips, etc. are made (yet!).
Can be used as a back-up for the ACD. 
<tr><td> Tkr1SSDVeto  
<td>F<td>   New version of the SSD Veto. For this variable, tracks which pass close
to a dead plane, buffer-saturated region, inter-wafer gap, gap between towers,
or dead strips, do not cause the veto count to be incremented if no cluster is found.
This almost certainly overdoes it: a more correct calculation would include
the probability for the track to cross the inactive region. (Coming soon!)
<tr><td> TkrVetoPlaneCrossed
<td>I<td>   Number of planes contributiong to the SSD Veto. This doesn't count the points
where a track crosses in a gap.
<tr><td> Tkr1CoreHC
<td>F<td>   Number of clusters within a roughly cylindrical region )(default radius 10 mm) 
around the TrackHits in each plane between the first and last on the best
track, excluding the clusters that belong to the track itself
<tr><td> TkrDispersion
<td>F<td>   the RMS of the distances between the 1st track and all others in the event.
For tracks which start "above" the first track distance is the 3-D distance
For the rest, distance is the doca of the head of the track to the axis of
the first track.
<tr><td> TkrUpstreamHC
<td>F<td>   The number of clusters in a cylinder (default radius 150 mm) up to 4 layers thick
above the head of the first track.
<tr><td> Tkr1CORERatio
<td>F<td>   the ratio of Tkr1CoreHC and Tkr1Hits
<tr><td> Tkr1LATEdge
<td>F<td>   Minimum distance to any LAT edge of the head of the best track
<tr><td> Tkr1FirstHits  
<td>F<td>   Number of initial TrackHits used to determine the starting direction  
<tr><td> Tkr1DifHits  
<td>F<td>   Difference between the number of x and y clusters associated with track  
<tr><td> Tkr1GhostFrac
<td>F<td>   Fraction of good hits on track 1 that are ghosts
<tr><td> Tkr1ToT255Frac
<td>F<td>   Fraction of good hits on track 1 that have ToT==255
<tr><td> Tkr1SaturatedFrac
<td>F<td>   Fraction of good hits on track 1 that are saturated
<tr><td> Tkr1WideFrac
<td>F<td>   Fraction of good hits on track 1 whose with is >= m_minWide (default=5)
<tr><td> Tkr1WiderFrac
<td>F<td>   Fraction of good hits on track 1 whose with is >= m_minWider (default=9)
<tr><td> Tkr1Gaps  
<td>F<td>   Total number of gaps in track  
<tr><td> Tkr1FirstGaps  
<td>F<td>   Number of gaps in first Tkr1FirstHits layers on track  
<tr><td> Tkr1DieEdge  
<td>F<td>   Distance from die (wafer) edge of initial point 
(0 is halfway between the dies, increases toward center of die)  
<tr><td> Tkr1KalThetaMS  
<td>F<td>   Multiple scattering angle (radians) referenced to first layer. 
The contributions from all the layers in the track are adjusted 
for the predicted energy in each layer, and weighted accordingly. 
So the result is sensitive to the particle type and 
the chosen energy-loss mechanism.       
<tr><td> Tkr1[X/Y/Z]Dir  
<td>D<td>   Track [x/y/z] direction cosine  
<tr><td> Tkr1Phi  
<td>F<td>   Track phi, radians 
(direction from which particle comes, not particle direction!) 
range: (0, 2pi)  
<tr><td> Tkr1Theta  
<td>F<td>   Track theta, radians (direction ditto)  
<tr><td> Tkr1[X/Y/Z]0  
<td>F<td>   Track [x/y/z] position at first hit  
</table>

@subsection ssdveto Diagnostic SSD Veto Variables (Optional, absent by default) 
@verbatim
(Turn on with jO parameter: ToolSvc.TkrValsTool.enableVetoDiagnostics = true;)
@endverbatim
<table>
<tr><th> Variable <th> Type  <th> Description                                
<tr><td> TkrVetoTrials
<td>I<td>   Difference between the plane number of the last plane crossed and the first, plus one
Any gaps above the last plane are not counted. (This may change soon.)
<tr><td> TkrVetoHitFound
<td>I<td>   Number of clusters found
<tr><td> TkrVetoUnknown
<td>I<td>   Missing clusters not ascribable to any inactive area
<tr><td> TkrVetoTower
<td>I<td>   Number of tower crossings (zero for now, because the propagator doesn't report them)
<tr><td> TkrVetoGapCorner
<td>I<td>   Number of missing clusters close to a wafer corner
<tr><td> TkrVetoGapEdge
<td>I<td>   Number of missing clusters close to a wafer edge
<tr><td> TkrVetoBadCluster
<td>I<td>   Number of missing clusters close to a dead strip
<tr><td> TkrVetoDeadPlane
<td>I<td>   Number of missing clusters close to a dead plane
<tr><td> TkrVetoTruncated
<td>I<td>   Number of missing clusters close to a truncated region
<tr><td> TkrVSSDVeto
<td>F<td>   Same as Tkr1SSDVeto, but for the veto track 
<tr><td> TkrVChisq
<td>F<td>   Same as Tkr1Chisq, but for the veto track
<tr><td> TkrVHits
<td>F<td>   Same as Tkr1Hits, but for the veto track
<tr><td> TkrVFirstLayer
<td>F<td>   Same as Tkr1FirstLayer, but for the veto track
<tr><td> TkrVKalEne
<td>F<td>   Same as Tkr1KalEne, but for the veto track
<tr><td> TkrVConEne
<td>F<td>   Same as Tkr1ConEne, but for the veto track
</table>

*/

StatusCode TkrCRValsTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;

    MsgStream log(msgSvc(), name());

    log << MSG::INFO  << "#################" << endreq << "# ";
    log << (m_useNew ? "New " : "Old ");
    log << "version" << endreq << "#################" << endreq;

    if((ValBase::initialize()).isFailure()) return StatusCode::FAILURE;

    m_messageCount = 0;

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

		 //Locate and store a pointer to the data service
        IService* iService = 0;
        if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
        {
           return sc;
        }
        m_dataSvc = dynamic_cast<DataSvc*>(iService);

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
        m_pAcdTool = 0;
        sc = toolSvc->retrieveTool("AcdValsTool", m_pAcdTool);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "Unable to find tool: " "AcdValsTool" << endreq;
            return sc;
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

    addItem("TkrCRNumTracks",   &TkrCR_Num_Tracks);
    addItem("TkrCRTrackLength", &TkrCR_TrackLength);

    addItem("TkrCR1Chisq",      &TkrCR_1_Chisq);
    addItem("TkrCR1FirstChisq", &TkrCR_1_FirstChisq);
    addItem("TkrCR1Hits",       &TkrCR_1_Hits);
    addItem("TkrCR1FirstHits",  &TkrCR_1_FirstHits);
    addItem("TkrCR1FirstLayer", &TkrCR_1_FirstLayer);
    addItem("TkrCR1LastLayer",  &TkrCR_1_LastLayer);
    addItem("TkrCR1DifHits",    &TkrCR_1_DifHits);

    addItem("TkrCR1ToT255Frac", &TkrCR_1_ToT255Frac);
    addItem("TkrCR1BothFrac",   &TkrCR_1_BothFrac);
    addItem("TkrCR1GhostFrac",  &TkrCR_1_GhostFrac);
    addItem("TkrCR1DiagFrac",   &TkrCR_1_DiagFrac);
    addItem("TkrCR1SaturatedFrac",  &TkrCR_1_SaturatedFrac);
    addItem("TkrCR1WideFrac",   &TkrCR_1_WideFrac);
    addItem("TkrCR1WiderFrac",  &TkrCR_1_WiderFrac);

    addItem("TkrCR1Gaps",       &TkrCR_1_Gaps);
    addItem("TkrCR1FirstGapPlane",&TkrCR_1_FirstGapPlane);
    addItem("TkrCR1XGap",       &TkrCR_1_GapX);
    addItem("TkrCR1YGap",       &TkrCR_1_GapY);
    addItem("TkrCR1FirstGaps",  &TkrCR_1_FirstGaps);

    addItem("TkrCR1Qual",       &TkrCR_1_Qual);
    addItem("TkrCR1Type",       &TkrCR_1_Type);
    addItem("TkrCR1TwrEdge",    &TkrCR_1_TwrEdge);
    addItem("TkrCR1PrjTwrEdge", &TkrCR_1_PrjTwrEdge);
    addItem("TkrCR1DieEdge",    &TkrCR_1_DieEdge);
    addItem("TkrCR1TwrGap",     &TkrCR_1_TwrGap);

    addItem("TkrCR1KalEne",     &TkrCR_1_KalEne);
    addItem("TkrCR1ConEne",     &TkrCR_1_ConEne);
    addItem("TkrCR1KalThetaMS", &TkrCR_1_KalThetaMS);

    addItem("TkrCR1XDir",       &TkrCR_1_xdir);
    addItem("TkrCR1YDir",       &TkrCR_1_ydir);
    addItem("TkrCR1ZDir",       &TkrCR_1_zdir);
    addItem("TkrCR1Phi",        &TkrCR_1_Phi);
    addItem("TkrCR1Theta",      &TkrCR_1_Theta);
    addItem("TkrCR1X0",         &TkrCR_1_x0);
    addItem("TkrCR1Y0",         &TkrCR_1_y0);
    addItem("TkrCR1Z0",         &TkrCR_1_z0);

    addItem("TkrCR1ToTFirst",   &TkrCR_1_ToTFirst);
    addItem("TkrCR1ToTAve",     &TkrCR_1_ToTAve);
    addItem("TkrCR1ToTTrAve",   &TkrCR_1_ToTTrAve);
    addItem("TkrCR1ToTAsym",    &TkrCR_1_ToTAsym);
    addItem("TkrCR1ChisqAsym",  &TkrCR_1_ChisqAsym);
    addItem("TkrCR1SSDVeto",    &TkrCR_1_SSDVeto, true);
    addItem("TkrCRPlaneCrossed",  &TkrCR_1_VetoPlaneCrossed);


    addItem("TkrCR1CoreHC",     &TkrCR_1_CoreHC);
    addItem("TkrCR1LATEdge",    &TkrCR_1_LATEdge);

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

StatusCode TkrCRValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    //Tkr_float = 5.5;
    //Tkr_int   = 123;

    //offset comes from Geometry
    double z0 = m_tkrGeom->gettkrZBot();

    //special stuff here
    TkrCR_1_FirstGapPlane = -1;
    
    //set TkrQueryClustersTool to return only normal clusters for the tests
    //  (even though the CR track may be a ghost!)
    pQueryClusters->setFilter(ITkrQueryClustersTool::NORMAL);

    double radThin  = m_tkrGeom->getAveConv(STANDARD); 
    double radThick = m_tkrGeom->getAveConv(SUPER); 
    double radTray  = m_tkrGeom->getAveRest(ALL);

    //Recover EventHeader Pointer
    //SmartDataPtr<Event::EventHeader> pEvent(m_pEventSvc, EventModel::EventHeader);

    // replace code below with the vector of pointers to all tracks  LSR
    // Recover Track associated info. 
    // NOTE: If no tracks are found ALL TKR variables are zero!  
    //SmartDataPtr<Event::TkrTrackCol>   
    //    pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol);

    // assemble the list CR Tracks

	  SmartDataPtr<Event::TkrTrackCol> 
        cRTrackCol(m_dataSvc, EventModel::TkrRecon::TkrCRTrackCol);

    if(!cRTrackCol) return sc;

    SmartDataPtr<Event::TkrClusterCol> 
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);


    // all variable values are preset to zero. 
    // Be sure to re-initialize the ones you care about  

    double die_width = m_tkrGeom->ladderPitch();
    int nDies = m_tkrGeom->nWaferAcross();

    if (cRTrackCol){   
        // Count number of tracks
//        int nTracks = pTracks->size();    RJ: doesn't work any more, with CR tracks in the mix

        int nTracks= cRTrackCol->size();
		TkrCR_Num_Tracks   = nTracks;
        if (cRTrackCol->empty()) return sc;
    
        // Get the first Track - it should be the "Best Track"
        Event::TkrTrackColConPtr pTrack = cRTrackCol->begin();
        const Event::TkrTrack* track_1 = *pTrack;

        TkrCR_1_Chisq        = track_1->getChiSquareSmooth();
        TkrCR_1_FirstChisq   = track_1->chiSquareSegment();
        TkrCR_1_FirstGaps    = track_1->getNumXFirstGaps() + track_1->getNumYFirstGaps();
        TkrCR_1_Qual         = track_1->getQuality();
        TkrCR_1_Type         = track_1->getStatusBits();
        TkrCR_1_Hits         = track_1->getNumFitHits();
        TkrCR_1_FirstHits    = track_1->getNumSegmentPoints();
        TkrCR_1_FirstLayer   = m_tkrGeom->getLayer(track_1->front()->getTkrId());
        TkrCR_1_LastLayer    = m_tkrGeom->getLayer(track_1->back()->getTkrId());
        TkrCR_1_Gaps         = track_1->getNumGaps();
        TkrCR_1_KalEne       = track_1->getKalEnergy(); 
        TkrCR_1_ConEne       = track_1->getInitialEnergy(); 
        TkrCR_1_KalThetaMS   = track_1->getKalThetaMS(); 
        TkrCR_1_RangeEne     = track_1->getRangeEnergy();
        TkrCR_1_DifHits      = track_1->getNumXHits()-track_1->getNumYHits();

        Point  x1 = track_1->getInitialPosition();
        Vector t1 = track_1->getInitialDirection();

        TkrCR_1_xdir        = t1.x();
        TkrCR_1_ydir        = t1.y();
        TkrCR_1_zdir        = t1.z();

        TkrCR_1_x0          = x1.x();
        TkrCR_1_y0          = x1.y();
        TkrCR_1_z0          = x1.z();

        // theta and phi are of direction of source, hence the minus sign
        // this code replaces atan and acos used before
        TkrCR_1_Phi         = (-t1).phi();
        if (TkrCR_1_Phi<0.0f) TkrCR_1_Phi += static_cast<float>(2*M_PI);
        TkrCR_1_Theta       = (-t1).theta();
        TkrCR_TrackLength = -(TkrCR_1_z0-z0)/TkrCR_1_zdir;

        double z_dist    = fabs((m_tkrGeom->trayHeight()+3.)/t1.z()); 
        double x_twr = globalToLocal(x1.x(), m_towerPitch, m_xNum);
        double y_twr = globalToLocal(x1.y(), m_towerPitch, m_yNum);

        double x_prj = x_twr - t1.x()*z_dist;
        double y_prj = y_twr - t1.y()*z_dist; 

        TkrCR_1_TwrEdge    = (fabs(x_twr) > fabs(y_twr)) ? fabs(x_twr) : fabs(y_twr);
        TkrCR_1_PrjTwrEdge = (fabs(x_prj) > fabs(y_prj)) ? fabs(x_prj) : fabs(y_prj);
        TkrCR_1_TwrEdge    = m_towerPitch/2. - TkrCR_1_TwrEdge;
        TkrCR_1_PrjTwrEdge = m_towerPitch/2. - TkrCR_1_PrjTwrEdge;

        // New section go compute gap lengths in tracker and cal
        double x_slope   = (fabs(t1.x()) > .0001)? t1.x():.00001;
        double s_x       = (sign(t1.x())*m_towerPitch/2. - x_twr)/x_slope; 
        double y_slope   = (fabs(t1.y()) > .0001)? t1.y():.00001;
        double s_y       = (sign(t1.y())*m_towerPitch/2. - y_twr)/y_slope;

        TkrCR_1_TwrGap = 0.; 
        if(s_x < s_y) { // Goes out x side of CAL Module
            if(x1.z() - z0 + s_x*t1.z() > 0) {
                double s_max = fabs((x1.z()-z0)/t1.z());
                TkrCR_1_TwrGap = gap/fabs(x_slope);
                if((TkrCR_1_TwrGap + s_x)> s_max ) TkrCR_1_TwrGap = s_max-s_x;
            }
        }
        else {          // Goes out y side
            if(x1.z() - z0 + s_y*t1.z() > 0) {
                double s_max = fabs((x1.z()-z0)/t1.z());
                TkrCR_1_TwrGap = gap/fabs(y_slope);
                if((TkrCR_1_TwrGap + s_y)> s_max ) TkrCR_1_TwrGap = s_max-s_y;
            }
        }

        // SSD Die loaction and edge... 
        double x_die = globalToLocal(x_twr, die_width, nDies);
        double y_die = globalToLocal(y_twr, die_width, nDies);


        TkrCR_1_DieEdge  = (fabs(x_die) > fabs(y_die)) ? fabs(x_die) : fabs(y_die);
        TkrCR_1_DieEdge  = die_width/2. - TkrCR_1_DieEdge; 

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

        // count up different types of hits for track 1 (same below for track 2)
        // count the number of real hits... may not be necessary but I'm nervous!
        int clustersOnTrack = 0;
        int toT255Count    = 0;
        int bothCount      = 0;
        int ghostCount     = 0;
        int diagCount      = 0;
        int saturatedCount = 0;
        int wideCount      = 0;
        int widerCount     = 0;
        int nToTs = 0; // count the good ToTs at the same time

        double mips;
        while(pHit != track_1->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            unsigned int bits = hit->getStatusBits();
            if((bits & Event::TkrTrackHit::HITISSSD)==0) continue;
            const Event::TkrCluster* cluster = hit->getClusterPtr();

            // maskZAPGHOSTS has all the ghost bits set, isSet() is true if any bit is set
            bool isMarked    = cluster->isSet(Event::TkrCluster::maskZAPGHOSTS);
            if(isMarked) {
                bool is255 = cluster->isSet(Event::TkrCluster::mask255);
                bool isGhost = cluster->isSet(Event::TkrCluster::maskGHOST);
                bool isDiagGhost = cluster->isSet(Event::TkrCluster::maskDIAGNOSTIC);
                if     (is255)                toT255Count++;
                else if(isGhost&&isDiagGhost) bothCount++;
                else if(isGhost)              ghostCount++;
                else if(isDiagGhost)          diagCount++;
            }
            mips = cluster->getMips();
            if(mips>-0.5) nToTs++;
            int rawToT = (int)cluster->ToT();
            if(rawToT==250) { saturatedCount++;}
            int width = (int)cluster->size();
            if(width>=m_minWide) wideCount++;
            if(width>=m_minWider) widerCount++;
            
            clustersOnTrack++;
        }
        if(clustersOnTrack>0) {
            TkrCR_1_ToT255Frac    = ((float)toT255Count)/clustersOnTrack;
            TkrCR_1_BothFrac      = ((float)bothCount)/clustersOnTrack;
            TkrCR_1_DiagFrac      = ((float)diagCount)/clustersOnTrack;
            TkrCR_1_GhostFrac     = ((float)ghostCount)/clustersOnTrack;
            TkrCR_1_SaturatedFrac = ((float)saturatedCount)/clustersOnTrack;
            TkrCR_1_WideFrac      = ((float)wideCount)/clustersOnTrack;
            TkrCR_1_WiderFrac     = ((float)widerCount)/clustersOnTrack;
        }

        // commented statements are part of the debug
        //std::vector<float> mipsVec(trackSize, 0.0);
        //std::vector<float> mipsBefore(trackSize, 0.0);
        //int seq = -1;
        //float maxMips = -1000;

        pHit = track_1->begin();
        //const Event::TkrTrackParams 
        //    params((*pHit)->getTrackParams(Event::TkrTrackHit::SMOOTHED));
        while(pHit != track_1->end()) {
            
            const Event::TkrTrackHit* hit = *pHit++;
            //seq++;
            unsigned int bits = hit->getStatusBits();

            int layer = m_tkrGeom->getLayer(hit->getTkrId());
            convType type = m_tkrGeom->getLayerType(layer);
            if (type==STANDARD)   {tkrTrackEnergy1 += cfThin;}
            else if (type==SUPER) {tkrTrackEnergy1 += cfThick;}

            // check if hit is in an ssd
            if ( !gapFound && (bits & Event::TkrTrackHit::HITISSSD)==0) {
                Point  gapPos = hit->getPoint(Event::TkrTrackHit::PREDICTED);
                TkrCR_1_GapX = gapPos.x();
                TkrCR_1_GapY = gapPos.y();
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
            double slope1  = fabs(hit->getMeasuredSlope(Event::TkrTrackHit::SMOOTHED));
            double slope = fabs(hit->getNonMeasuredSlope(Event::TkrTrackHit::SMOOTHED));

            // theta1 is the projected angle across the strip
            double theta1       = atan(slope1);

            double aspectRatio = 0.228/0.400;
            double threshold   =  0.25;   // Mips
            //double countThreshold = 15.0; // counts
            //double normFactor  =  1./53.;

            double mips = cluster->getMips();
            //mipsBefore[seq] = mips;
            if(mips<-0.5) continue;

            double path1 = 1.0;

            // get the path length for the hit
            // tries to get the average length
            // the calculation is part analytic, part approximation and part fudge.
            //   more work is definitely in order!

            // theta1 first
            double costh1 = cos(theta1);
            if (size==1) {
                double sinth1 = sin(theta1);
                if (slope1< aspectRatio) {
                    path1 = (1./costh1*(aspectRatio-slope1) + 
                        (1/costh1 - 0.5*threshold)*(2*threshold*sinth1))
                        /(aspectRatio - slope1 + 2*threshold*sinth1);
                } else if (slope1<aspectRatio/(1-2.*threshold*costh1)) {
                    path1 = 1.0; //1/costh1 - threshold*costh1;
                } else { 
                    path1 = 1.0;
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
            //mipsVec[seq] = mips;

            if(mips > max_ToT) max_ToT = mips; 
            if(mips < min_ToT) min_ToT = mips; 
            hit_counter++;  
            if (hit_counter==1) TkrCR_1_ToTFirst = mips;
            TkrCR_1_ToTAve += mips;
            // first 2 valid clusters
            if(hit_counter < 3) {
                first_ToTs += mips;
                chisq_first += hit->getChiSquareSmooth();
            }
            // last 2 valid clusters
            if(hit_counter > nToTs-2){
                last_ToTs += mips;
                chisq_last += hit->getChiSquareSmooth();
            }
        }

        // need 4 clusters to calculate an asymmetry
        if(nToTs>3&&first_ToTs+last_ToTs>0) {
                TkrCR_1_ToTAsym = (last_ToTs - first_ToTs)/(first_ToTs + last_ToTs);
        }

        TkrCR_1_ToTAve /= std::max(1, nToTs);

        // at least three tracks to do a real truncated mean, if not use normal mean
        if (nToTs>2) {
            TkrCR_1_ToTTrAve = (TkrCR_1_ToTAve*nToTs - max_ToT - min_ToT)/(nToTs-2.);
        } else { 
            TkrCR_1_ToTTrAve = TkrCR_1_ToTAve;
        }

        tkrTrackEnergy1 /= fabs(TkrCR_1_zdir);


        TkrCR_1_FirstGapPlane = gapId; 

        // Chisq Asymmetry - Front vs Back ends of tracks
        if (chisq_last+chisq_first>0) TkrCR_1_ChisqAsym = 
            (chisq_last - chisq_first)/(chisq_last + chisq_first);


        int firstPlane = m_tkrGeom->getPlane(track_1->front()->getTkrId()); 
        int firstLayer = m_tkrGeom->getLayer(firstPlane);
        double zFirstLayer = m_tkrGeom->getLayerZ(firstLayer);

        double costh = fabs(t1.z());
        double secth = 1./costh;

        // for the footprints
        double secthX = 1./sqrt(1.0 - t1.x()*t1.x());
        double secthY = 1./sqrt(1.0 - t1.y()*t1.y());

        //double xVetoRgn = _vetoRegion*secthX;
        //double yVetoRgn = _vetoRegion*secthY;

        // SSD Veto stuff here:
        // First Track stuff as before
        TkrCR_1_SSDVeto = SSDEvaluation(track_1); 

        TkrCR_1_VetoPlaneCrossed = (int)floor(m_VetoPlaneCrossed + 0.5); 
        TkrCR_1_VetoTrials       = (int)floor(m_VetoTrials + 0.5);
        TkrCR_1_VetoUnknown      = (int)floor(m_VetoUnknown + 0.5);
        TkrCR_1_VetoDeadPlane    = (int)floor(m_VetoDeadPlane + 0.5);
        TkrCR_1_VetoTruncated    = (int)floor(m_VetoTruncated + 0.5); 
        TkrCR_1_VetoTower        = (int)floor(m_VetoTower + 0.5);   
        TkrCR_1_VetoGapCorner    = (int)floor(m_VetoGapCorner + 0.5);
        TkrCR_1_VetoGapEdge      = (int)floor(m_VetoGapEdge + 0.5);   
        TkrCR_1_VetoBadCluster   = (int)floor(m_VetoBadCluster + 0.5);

        // minimum distance from any edge, measured from the edge of the active area
        double deltaEdge = 0.5*(m_towerPitch - m_tkrGeom->trayWidth()) 
            - m_tkrGeom->siDeadDistance() ;
        double tkrXLo = m_tkrGeom->getLATLimit(0,LOW)  + deltaEdge;
        double tkrXHi = m_tkrGeom->getLATLimit(0,HIGH) - deltaEdge;
        double tkrYLo = m_tkrGeom->getLATLimit(1,LOW)  + deltaEdge;
        double tkrYHi = m_tkrGeom->getLATLimit(1,HIGH) - deltaEdge;

        double xEdge = std::min(TkrCR_1_x0-tkrXLo, tkrXHi- TkrCR_1_x0);
        double yEdge = std::min(TkrCR_1_y0-tkrYLo, tkrYHi- TkrCR_1_y0);

        TkrCR_1_LATEdge = (float) std::min(xEdge, yEdge);

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
            TkrCR_1_CoreHC += coreHits;
        }
	}
	return sc;
}

double TkrCRValsTool::towerEdge(Point pos) const
{
    double edge = 0.; 
    double x = pos.x();
    double y = pos.y();

    double x_twr = globalToLocal(x, m_towerPitch, m_xNum);
    double y_twr = globalToLocal(y, m_towerPitch, m_yNum);

    edge = 0.5*m_towerPitch - std::max(fabs(x_twr),fabs(y_twr));
    return edge;
}

float TkrCRValsTool::SSDEvaluation(const Event::TkrTrack* track) 
{// Method to compute the number of SSD Vetoes for the given track

    MsgStream log(msgSvc(), name());

    Point  x1 = track->getInitialPosition();
    Vector t1 = track->getInitialDirection();

    Event::TkrTrackHitVecConItr pHit = track->begin();
    const Event::TkrTrackParams params((*pHit)->getTrackParams(Event::TkrTrackHit::SMOOTHED));


    float Tkr_ConEne       = track->getInitialEnergy(); 

    int topPlane = m_tkrGeom->numPlanes()-1; 
    double topOfTkr = m_tkrGeom->getPlaneZ(topPlane) + 1.0;
    double arc_min = fabs((topOfTkr-x1.z())/t1.z());
    arc_min = std::min( arc_min, maxPath);

    bool upwards = true;

    m_VetoPlaneCrossed = _badFloat; 
    m_VetoTrials = _badFloat;
    m_SSDVeto = _badFloat;
    m_VetoUnknown = _badFloat;
    m_VetoDeadPlane = _badFloat;
    m_VetoTruncated = _badFloat; 
    m_VetoTower = -1.;   
    m_VetoGapCorner = _badFloat;
    m_VetoGapEdge = _badFloat;   
    m_VetoBadCluster = _badFloat;
    m_VetoHitFound = _badFloat;

    try {
        m_G4PropTool->setStepStart(params, x1.z(), upwards);
        m_G4PropTool->step(arc_min);
    } catch( std::exception& /*e*/) {
        printHeader(log);
        setAnaTupBit();
        log << "See previous exception printout." << endreq;
        log << " Skipping the TKR Veto calculations" << endreq;
        return m_SSDVeto;

    } catch (...) {
        printHeader(log);
        setAnaTupBit();
        log << "Unknown exception, see previous exception message, if any" << endreq;
        log << "Skipping the TKR Veto calculations" << endreq;
        log << "Initial track parameters: pos: " << x1 << endreq 
            << "dir: " << t1 << " arclen: " << arc_min << endreq;
        return m_SSDVeto;
    }

    m_VetoPlaneCrossed = 0.0; 
    m_VetoTrials = 0.0;
    m_SSDVeto = 0.0;
    m_VetoUnknown = 0.0;
    m_VetoDeadPlane = 0.0;
    m_VetoTruncated = 0.0; 
    m_VetoTower = -1.;   
    m_VetoGapCorner = 0.0;
    m_VetoGapEdge = 0.0;   
    m_VetoBadCluster = 0.0;
    m_VetoHitFound = 0.0;

    int numSteps = m_G4PropTool->getNumberSteps();

    idents::VolumeIdentifier volId;
    idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();

    //int firstPlane = m_tkrGeom->getPlane(track->front()->getTkrId()); 
    //int firstLayer = m_tkrGeom->getLayer(firstPlane);
    //double zFirstLayer = m_tkrGeom->getLayerZ(firstLayer);

    //double costh = fabs(t1.z());
    //double secth = 1./costh;

    // for the footprints
    double secthX = 1./sqrt(1.0 - t1.x()*t1.x());
    double secthY = 1./sqrt(1.0 - t1.y()*t1.y());

    double xVetoRgn = _vetoRegion*secthX;
    double yVetoRgn = _vetoRegion*secthY;

    // Note: skip the first vol - as this is the head SSD

    double arcLen = m_G4PropTool->getStepArcLen(0);
    bool isFirstPlane = true;
    int previousplane = -1;
    for(int istep = 1; istep < numSteps; ++istep) { 
        volId = m_G4PropTool->getStepVolumeId(istep);
        volId.prepend(prefix);
        Point x_step       = m_G4PropTool->getStepPosition(istep);
        arcLen += m_G4PropTool->getStepArcLen(istep);

        bool forward = false;
        const Event::TkrTrackParams newParams = 
            m_G4PropTool->getTrackParams(arcLen, Tkr_ConEne, forward);

        //std::cout << "pos " << x_step << ", step " << istep << ", arcLen " << arcLen << std::endl;
        //std::cout << "err " << sqrt(newParams(xPosIdx,xPosIdx)) << " "
        //        << sqrt(newParams(yPosIdx,yPosIdx)) << std::endl;
        // we're outside the LAT
        if(x_step.z() > topOfTkr || !m_tkrGeom->isInActiveLAT(x_step) ) break; 

        // check that it's really a TKR hit (probably overkill)
        if(volId.size() != 9) continue; 
        if(!(volId[0]==0 && volId[3]==1)) continue; // !(LAT && TKR)
        if(volId[6]> 1) continue;  //It's a converter or some other tray element!
 
        // now check if there's a hit near the extrapolated track!
        // if there is, reset the veto counter... we want leading non-hits
        //std::cout << "Tkr1SSDVeto volId " << volId.name() << std::endl;

        //int tower = idents::TowerId(volId[2], volId[1]).id();
        int tray = volId[4];
        int view = volId[5];
        int face = volId[6];
        int layer = m_tkrGeom->trayToBiLayer(tray, face);
        int plane = m_tkrGeom->trayToPlane(tray, face);
        if(previousplane==plane) continue; // already looked at this plane
        previousplane=plane;
        m_VetoPlaneCrossed++; //

        int firstPlane = -1;
        if(isFirstPlane) {
            isFirstPlane = false;
            firstPlane = plane;
        }

        // I think we want to do this, most likely this is missed for
        //   some good reason.
        // on the other hand, it is a missed hit, so remove test for the new code


        m_VetoTrials = abs(plane-firstPlane) + 1;

        //if(!m_useNew&&layer==firstLayer) continue;

        double vetoRgn = (view==0 ? xVetoRgn : yVetoRgn);

        int nVetoHits = pQueryClusters->numberOfHitsNear(view, layer, 
            vetoRgn, x_step, t1);

        if (nVetoHits>0) {
            // found a hit, reset the SSDVeto
            m_SSDVeto = 0.0;
            if(m_enableVetoDiagnostics) { m_VetoHitFound++; }

        } else { 
            // no hit
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
                m_SSDVeto += 1.0; 
            }

            if(m_enableVetoDiagnostics) {
                switch (stage) 
                {
                case -1:
                    m_VetoUnknown++;
                    break;
                case 1: 
                    m_VetoDeadPlane++;
                    break;
                case 2:
                    m_VetoTruncated++;
                    break;
                case 3:
                    m_VetoTower++;
                    break;
                case 4:
                    m_VetoGapCorner++;
                    break;
                case 5:
                    m_VetoGapEdge++;
                    break;
                case 6:
                    m_VetoBadCluster++;
                    break;
                default:
                    std::cout << "shouldn't get here, stage = " 
                        << stage << std::endl;
                } // end switch
            } // end diagnostics
        } // end no hit 
    } // end loop over steps

    return m_SSDVeto;
}
