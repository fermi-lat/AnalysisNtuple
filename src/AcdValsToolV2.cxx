
/** @file AcdValsTool.cxx
@brief Calculates the Adc analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/AcdValsToolV2.cxx,v 1.3 2009/12/21 20:07:27 jrb Exp $
*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"   
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "Event/Recon/AcdRecon/AcdReconV2.h"
#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"
#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdEventTopology.h"


#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Digi/AcdDigi.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "AcdUtil/AcdTileFuncs.h"
//#include "AcdRecon/AcdGap.h" (Don't include AcdRecon, use ints instead of enums)

#include <algorithm>
#include <numeric>
#include <map>

/** @class Acd2ValsTool
@brief Calculates Acd Values
*/

class Acd2ValsTool : public ValBase

{

public:

  Acd2ValsTool( const std::string& type, 
		const std::string& name, 
		const IInterface* parent);
  
  virtual ~Acd2ValsTool() { }
  
  StatusCode initialize();
  
  StatusCode calculate();
  
protected:

  void loadTopologyVars(const Event::AcdReconV2& pACD);

  void reconId(const Event::AcdReconV2* pACD);
  
  void setId(const std::vector<Event::AcdTkrHitPoca*> trackUp,
	     unsigned int &retId, bool findRibbon=false);
  
  void findId(const std::vector<Event::AcdTkrHitPoca*>& vec, 
	      idents::AcdId &retId, bool findRibbon=false);  
  
private:
  
  // Global ACD Tuple Items
  unsigned int ACD_Tile_Count; 
  unsigned int ACD_Ribbon_Count;
  unsigned int ACD_Veto_Count; 
  unsigned int ACD_Veto_Faces;
  float        ACD_Total_Energy;
  float        ACD_Total_Ribbon_Energy;

  unsigned int ACD_TileIdRecon;
  unsigned int ACD_RibbonIdRecon;
  int          ACD_ActiveDist_TrackNum;
  
  // Variables computed by looping over all tracks w.r.t. hit tiles
  float ACD_ActiveDist3D;
  float ACD_ActiveDist3D_Err;
  float ACD_ActiveDist_Energy;  
  
  // Variables computed by looping over all tracks w.r.t. hit ribbons
  float ACD_ribbon_ActiveDist;
  float ACD_ribbon_ActiveDist_Err;
  float ACD_ribbon_ActiveLength;
  float ACD_ribbon_EnergyPmtA;
  float ACD_ribbon_EnergyPmtB;
  
  // Variables computed by looping over all tracks w.r.t. gaps in the ACD
  float ACD_Corner_DOCA;
  float ACD_TkrHole_Dist;
  float ACD_TkrRibbon_Dist; 
  float ACD_TkrRibbon_Dist_Err; 
  float ACD_TkrRibbonLength; 
  
  // Variables computed by taking best track w.r.t. hit tiles
  float ACD_Tkr1ActiveDist;
  float ACD_Tkr1ActiveDist_Err;
  float ACD_Tkr1ActiveDist_Energy;
  float ACD_Tkr1ActiveDist_Arc;
  float ACD_Tkr1Energy15;
  float ACD_Tkr1Energy30;
  float ACD_Tkr1Energy45;
  
  // Variables computed by taking best track w.r.t. hit ribbons
  float ACD_Tkr1_ribbon_ActiveDist;
  float ACD_Tkr1_ribbon_ActiveDist_Err;
  float ACD_Tkr1_ribbon_ActiveLength;
  float ACD_Tkr1_ribbon_EnergyPmtA;
  float ACD_Tkr1_ribbon_EnergyPmtB;
  
  // Variables computed by taking best w.r.t. gaps in the ACD    
  float ACD_Tkr1Corner_DOCA;
  float ACD_Tkr1Hole_Dist;
  float ACD_Tkr1Ribbon_Dist;
  float ACD_Tkr1Ribbon_Dist_Err;
  float ACD_Tkr1RibbonLength;
    
  // New stuff EAC
  float ACD_Tkr1_VetoSigmaHit;
  float ACD_Tkr1_VetoSigmaGap;
  float ACD_Tkr_VetoSigmaHit;
  float ACD_Tkr_VetoSigmaGap;
  
  // Variables about number of ACD tiles by row
  unsigned int ACD_tileTopCount;
  unsigned int ACD_tileCount0;
  unsigned int ACD_tileCount1;
  unsigned int ACD_tileCount2;
  unsigned int ACD_tileCount3;
  
  unsigned int ACD_countRow3Readout;
  
  float ACD_energyTop;
  float ACD_energyRow0;
  float ACD_energyRow1;
  float ACD_energyRow2;
  float ACD_energyRow3;
  
  // Services
  IGlastDetSvc *m_detSvc;
  
  // Algorithm parameters
  double m_vetoThresholdMeV;

  // Prefix to add to columns in tuple [Acd] (production) or [Acd2] (regression testing)
  std::string m_prefix;

};

/** @page anatup_vars
@section acd2valstool AcdValsV2Tool Variables
Notes
- Default Doca/ActiveDistance is -2000.
- Active distance is negative if a track is outside a tile, 
positive if inside.

<table>
<tr><th> Variable <th> Type <th> Description					
<tr><td> Acd2TileCount 
<td>U<td>   Number of tiles fired
<tr><td> Acd2RibbonCount 
<td>U<td>   Number of ribbons fired
<tr><td> Acd2VetoCount 
<td>U<td>   Number of vetos fired
<tr><td> Acd2VetoFaces
<td>U<td>   Number of Acd faces with vetos
<tr><td> Acd2TotalEnergy 	
<td>F<td>   Total energy deposited in ACD Tiles
<tr><td> Acd2RibbonEnergy 	
<td>F<td>   Total energy deposited in ACD Ribbons

<tr><td> Acd2TileIdRecon
<td>U<td> Tile identifier that was pierced by the reconstructed track.  
A value of 899 (N/A) is the default and denotes that no ACD tile was 
intersected by a reconstructed track.
<tr><td> Acd2RibbonIdRecon
<td>U<td> Ribbon identifier that was pierced by the reconstructed track.  
A value of 899 (N/A) is the default and denotes that no ACD ribbon was 
intersected by a reconstructed track.


<tr><td> Acd2ActiveDist3D  	
<td>F<td>   Active Distance most likely to give a veto.  
Corresponds to Act. Dist. that is greater than a set 
an energy dep. min. distance and has the largest pulse height
<tr><td> Acd2ActiveDist3DErr
<td>F<td>   Error on most likely veto active distance of any track 
to the edge of any tile 
<tr><td> Acd2ActDistTileEnergy 
<td>F<td>   The deposited energy in the corresponding hit tile 
<tr><td> Acd2ActDistTrackNum
<td>F<td>   Track number of track which was used for AcdActiveDist3D. 
Track numbering starts at zero; best track number is zero; -1 means no track

<tr><td> Acd2RibbonActDist
<td>F<td>   Largest active distance to any ribbon 
(considered as a straight at the center of the ribbon) 
<tr><td> Acd2RibbonActDistErr  
<td>F<td>   Error on the smallest active distance to any ribbon 
(considered as a straight line of no thickness) 
<tr><td> Acd2RibbonActLength
<td>F<td>   Length along ribbon where point of closest approach occured. 
0 is center of ribbon; + going towards +x or +y side of ACD
<tr><td> Acd2RibbonActEnergyPmtA
<td>F<td>   The deposited energy in the A PMT of the corresponding hit ribbon
<tr><td> Acd2RibbonActEnergyPmtB

<td>F<td>   The deposited energy in the B PMT of the corresponding hit ribbon
<tr><td> Acd2CornerDoca 
<td>F<td>   Minimum Distance of Closest Approach of best track to the corner side gaps
This variable is signed to match the direction of the overlaps in the tiles. 
The gap appears larger for tracks coming from the + side than the - side.  
<tr><td> Acd2TkrHoleDist (fixme)
<td>F<td>   Minimum Distance of Closest Approach of any track to any of the tile screw holes
<tr><td> Acd2TkrRibbonDist
<td>F<td>   Minimum Distance of Closest Approach of any track to any ribbons that cover gaps
<tr><td> Acd2TkrRibbonDistErr
<td>F<td>   Error On Minimum Distance of Closest Approach of any track to any ribbons that cover gaps
<tr><td> Acd2TkrRibbonLength
<td>F<td>   Length along ribbon where point of closest approach occured.
0 is center of ribbon + going towards +x or +y side of ACD.


<tr><td> Acd2Tkr1ActiveDist 
<td>F<td>   Largest active distance from  track 1 to the edge of any tile
<tr><td> Acd2Tkr1ActiveDistErr
<td>F<td>   Error on largest active distance from track 1 to the edge of any tile
<tr><td> Acd2Tkr1ActDistTileEnergy
<td>F<td>   The deposited energy in the corresponding hit tile
<tr><td> Acd2Tkr1ActDistTileArc
<td>F<td>   Length from head of track to point where active distances was calculated
<tr><td> Acd2Tkr1Energy15
<td>F<td>   Energy in 15 deg. cone ahead of Track
<tr><td> Acd2Tkr1Energy30
<td>F<td>   Energy in 45 deg. cone ahead of Track
<tr><td> Acd2Tkr1Energy45
<td>F<td>   Energy in 45 deg. cone ahead of Track

<tr><td> Acd2Tkr1RibbonActDist   
<td>F<td>   Largest active distance to any ribbon 
(considered as a straight line at the center of the ribbon) 
<tr><td> Acd2Tkr1RibbonActDistErr   
<td>F<td>   Error on the smallest active distance to any ribbon 
<tr><td> Acd2Tkr1RibbonActLength
<td>F<td>   Length along ribbon where point of closest approach occured. 
0 is center of ribbon + going towards +x or +y side of ACD
<tr><td> Acd2Tkr1RibbonActEnergyPmtA   
<td>F<td>   The deposited energy in the A PMT of the corresponding hit ribbon
<tr><td> Acd2Tkr1RibbonActEnergyPmtB   
<td>F<td>   The deposited energy in the A PMT of the corresponding hit ribbon

<tr><td> Acd2Tkr1CornerDoca 
<td>F<td>   Minimum Distance of Closest Approach of best track to the corner side gaps 
This variable is signed to match the direction of the overlaps in the tiles. 
The gap appears larger for tracks coming from the + side than the - side.  
<tr><td> Acd2Tkr1HoleDist (fixme)
<td>F<td>   Minimum Distance of Closest Approach to best track to any of the tile screw holes
<tr><td> Acd2Tkr1RibbonDist
<td>F<td>   Minimum Distance of Closest Approach to best track to any ribbons that cover gaps 
<tr><td> Acd2Tkr1RibbonDistErr
<td>F<td>   Error on Minimum Distance of Closest Approach to best track to any ribbons that cover gaps 
<tr><td> Acd2Tkr1RibbonLength 
<td>F<td>   Length along ribbon where point of closest approach occured. 
0 is center of ribbon; + going towards +x or +y side of ACD

<tr><td> Acd2Tkr1VetoSigmaHit
<td>F<td>Number of sigmas less than an expected mip for signal combined with the
number of sigmas track propagation is away from tile or ribbon most likely to veto the best track.
<tr><td> Acd2Tkr1VetoSigmaGap
<td>F<td>Number of sigmas track propagation is away from clostest GAP in ACD for the best track.
<tr><td> Acd2TkrVetoSigmaHit
<td>F<td>Number of sigmas less than an expected mip for signal combined with the
number of sigmas track propagation is away from tile or ribbon most likely to veto any track.
<tr><td> Acd2TkrVetoSigmaGap
<td>F<td>Number of sigmas track propagation is away from closest GAP in ACD for any track.

<tr><td> Acd2NoTop
<td>U<td>   Hit Tile counts top of ACD
<tr><td> Acd2NoSideRow[0...3] 
<td>U<td>   Hit Tile counts for side row [0...3] that have energy > TileCountThreshold (= .8)
<tr><td> Acd2NoRow3Readout
<td>U<td>   Hit Tile counds for side row 3, no threshold cut
<tr><td> Acd2EnergyTop
<td>F<td>   Total energy deposited in top tiles
<tr><td> Acd2EnergyRow[0...3]
<td>F<td>   Total energy deposited in the tiles in side row 0 to 3
</table>
*/


// Static factory for instantiation of algtool objects
static ToolFactory<Acd2ValsTool> s_factory;
const IToolFactory& Acd2ValsToolFactory = s_factory;

// Standard Constructor
Acd2ValsTool::Acd2ValsTool(const std::string& type, 
			   const std::string& name, 
			   const IInterface* parent)
  : ValBase( type, name, parent )
{    
  // Declare additional interface
  declareInterface<IValsTool>(this); 
  
  // Threshold in MeV
  declareProperty("VetoThresholdMeV", m_vetoThresholdMeV=0.0);

  // Prefix for tuple column names
  declareProperty("Prefix", m_prefix="Acd2");
}


StatusCode Acd2ValsTool::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;
  
  MsgStream log(msgSvc(), name());
  
  if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;
  
  setProperties();
  // get the services
  
  // use the pointer from ValBase  
  if( serviceLocator() ) {
    
    // find GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc, true).isFailure()){
      log << MSG::INFO << "Couldn't find the GlastDetSvc!" << endreq;
      log << MSG::INFO << "Will be unable to calculate ACD_TkrHitsCount" << endreq;
      m_vetoThresholdMeV = 0.4;
    } else {
      StatusCode sc = m_detSvc->getNumericConstByName("acd.vetoThreshold", 
						      &m_vetoThresholdMeV);
      if (sc.isFailure()) {
	log << MSG::INFO << 
	  "Unable to retrieve threshold, setting the value to 0.4 MeV" << endreq;
	m_vetoThresholdMeV = 0.4;
      }
    }
    
  } else {
    return StatusCode::FAILURE;
  }
  
  // load up the map
  addItem(m_prefix + "TileCount",    &ACD_Tile_Count);
  addItem(m_prefix + "RibbonCount", &ACD_Ribbon_Count);
  addItem(m_prefix + "VetoCount",    &ACD_Veto_Count);
  addItem(m_prefix + "VetoFaces", &ACD_Veto_Faces);
  addItem(m_prefix + "TotalEnergy", &ACD_Total_Energy);
  addItem(m_prefix + "RibbonEnergy", &ACD_Total_Ribbon_Energy);
  addItem(m_prefix + "TileIdRecon", &ACD_TileIdRecon);
  addItem(m_prefix + "RibbonIdRecon", &ACD_RibbonIdRecon);
  
  addItem(m_prefix + "ActiveDist3D",   &ACD_ActiveDist3D);
  addItem(m_prefix + "ActiveDist3DErr",   &ACD_ActiveDist3D_Err);
  addItem(m_prefix + "ActDistTileEnergy",   &ACD_ActiveDist_Energy);
  addItem(m_prefix + "ActDistTrackNum", &ACD_ActiveDist_TrackNum);
  
  addItem(m_prefix + "RibbonActDist", &ACD_ribbon_ActiveDist);
  addItem(m_prefix + "RibbonActDistErr", &ACD_ribbon_ActiveDist_Err);
  addItem(m_prefix + "RibbonActLength", &ACD_ribbon_ActiveLength);
  addItem(m_prefix + "RibbonActEnergyPmtA", &ACD_ribbon_EnergyPmtA);
  addItem(m_prefix + "RibbonActEnergyPmtB", &ACD_ribbon_EnergyPmtB);

  addItem(m_prefix + "CornerDoca",    &ACD_Corner_DOCA);
  addItem(m_prefix + "TkrHoleDist",    &ACD_TkrHole_Dist);
  addItem(m_prefix + "TkrRibbonDist",    &ACD_TkrRibbon_Dist);
  addItem(m_prefix + "TkrRibbonDistErr",    &ACD_TkrRibbon_Dist_Err);
  addItem(m_prefix + "TkrRibbonLength",    &ACD_TkrRibbonLength);
  
  addItem(m_prefix + "Tkr1ActiveDist", &ACD_Tkr1ActiveDist);
  addItem(m_prefix + "Tkr1ActiveDistErr", &ACD_Tkr1ActiveDist_Err);
  addItem(m_prefix + "Tkr1ActiveDistArc", &ACD_Tkr1ActiveDist_Arc);
  addItem(m_prefix + "Tkr1ActDistTileEnergy", &ACD_Tkr1ActiveDist_Energy);
  addItem(m_prefix + "Tkr1Energy15", &ACD_Tkr1Energy15);
  addItem(m_prefix + "Tkr1Energy30", &ACD_Tkr1Energy30);
  addItem(m_prefix + "Tkr1Energy45", &ACD_Tkr1Energy45);
  
  addItem(m_prefix + "Tkr1RibbonActDist", &ACD_Tkr1_ribbon_ActiveDist);
  addItem(m_prefix + "Tkr1RibbonActDistErr", &ACD_Tkr1_ribbon_ActiveDist_Err);
  addItem(m_prefix + "Tkr1RibbonActLength", &ACD_Tkr1_ribbon_ActiveLength);
  addItem(m_prefix + "Tkr1RibbonActEnergyPmtA", &ACD_Tkr1_ribbon_EnergyPmtA);
  addItem(m_prefix + "Tkr1RibbonActEnergyPmtB", &ACD_Tkr1_ribbon_EnergyPmtB);
  
  addItem(m_prefix + "Tkr1CornerDoca",    &ACD_Tkr1Corner_DOCA);
  addItem(m_prefix + "Tkr1HoleDist",    &ACD_Tkr1Hole_Dist);
  addItem(m_prefix + "Tkr1RibbonDist",    &ACD_Tkr1Ribbon_Dist);
  addItem(m_prefix + "Tkr1RibbonDistErr",    &ACD_Tkr1Ribbon_Dist_Err);
  addItem(m_prefix + "Tkr1RibbonLength",    &ACD_Tkr1RibbonLength);    
  
  addItem(m_prefix + "NoTop",        &ACD_tileTopCount);
  addItem(m_prefix + "NoSideRow0",   &ACD_tileCount0);
  addItem(m_prefix + "NoSideRow1",   &ACD_tileCount1);
  addItem(m_prefix + "NoSideRow2",   &ACD_tileCount2);   
  addItem(m_prefix + "NoSideRow3",   &ACD_tileCount3);
  
  addItem(m_prefix + "NoRow3Readout", &ACD_countRow3Readout);
  addItem(m_prefix + "EnergyTop",  & ACD_energyTop);
  addItem(m_prefix + "EnergyRow0", & ACD_energyRow0);
  addItem(m_prefix + "EnergyRow1", & ACD_energyRow1);
  addItem(m_prefix + "EnergyRow2", & ACD_energyRow2);
  addItem(m_prefix + "EnergyRow3", & ACD_energyRow3);
  
  addItem(m_prefix + "Tkr1VetoSigmaHit", & ACD_Tkr1_VetoSigmaHit);
  addItem(m_prefix + "Tkr1VetoSigmaGap", & ACD_Tkr1_VetoSigmaGap);
  addItem(m_prefix + "TkrVetoSigmaHit", & ACD_Tkr_VetoSigmaHit);
  addItem(m_prefix + "TkrVetoSigmaGap", & ACD_Tkr_VetoSigmaGap);   
  
  zeroVals();
  ACD_ActiveDist_TrackNum = -1;
  
  return sc;
}


StatusCode Acd2ValsTool::calculate()
{

  StatusCode sc = StatusCode::SUCCESS;
  
  // Recover pointers to ACD Recon results
  SmartDataPtr<Event::AcdReconV2> pACD(m_pEventSvc,EventModel::AcdReconV2::Event);

  // Recover pointers to CalClusters and Xtals(not currently used)
  //SmartDataPtr<Event::CalClusterCol>   pCals(m_pEventSvc,EventModel::CalRecon::CalClusterCol);

  // Recover Track associated info. (not currently used)
  //SmartDataPtr<Event::TkrFitTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
  //SmartDataPtr<Event::TkrVertexCol>    pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);

  //Make sure we have valid ACD data
  if (!pACD) {
    return StatusCode::FAILURE;
  }

  loadTopologyVars(*pACD); 

  reconId(pACD);

  // Make a map relating AcdId to calibrated hits
  std::map<idents::AcdId, const Event::AcdHit*> hitMap;
    
  const Event::AcdHitCol& hitCol = pACD->getHitCol();
  int nHit = hitCol.size();
  
  // Loop over the hits and fill the maps
  for (int iHit(0); iHit < nHit; iHit++ ){
    const Event::AcdHit* aHit = hitCol[iHit];
    const idents::AcdId& id = aHit->getAcdId();
    hitMap[id] = aHit;
  }

  static const float maxSigma = 10000.;
  static const float maxSigmaSq = maxSigma*maxSigma;
  static const float maxDist = 2000.;
  static const float negMaxDist = -2000.;
  
  // Reset variables for loop over all tracks
  ACD_ActiveDist3D = ACD_ribbon_ActiveDist = negMaxDist;
  ACD_ActiveDist_Energy = ACD_ribbon_EnergyPmtA = ACD_ribbon_EnergyPmtB = 0.;
  ACD_ActiveDist3D_Err = ACD_ribbon_ActiveDist_Err = ACD_TkrRibbon_Dist_Err = -1.;	
  ACD_Corner_DOCA = ACD_TkrRibbon_Dist = ACD_TkrHole_Dist = negMaxDist;
  ACD_ribbon_ActiveLength = ACD_TkrRibbonLength = -10000.;
  
  // Reset variables for best track
  ACD_Tkr1ActiveDist = ACD_Tkr1_ribbon_ActiveDist = negMaxDist;
  ACD_Tkr1ActiveDist_Energy = ACD_Tkr1_ribbon_EnergyPmtA 
    = ACD_Tkr1_ribbon_EnergyPmtB = 0.;
  ACD_Tkr1ActiveDist_Err = ACD_Tkr1_ribbon_ActiveDist_Err = ACD_Tkr1Ribbon_Dist_Err = -1.;
  ACD_Tkr1ActiveDist_Arc = -1.;
  ACD_Tkr1Corner_DOCA = ACD_Tkr1Ribbon_Dist = ACD_Tkr1Hole_Dist = negMaxDist;
  ACD_Tkr1_ribbon_ActiveLength = ACD_Tkr1RibbonLength = -10000.;
  ACD_Tkr1Energy15 = ACD_Tkr1Energy30 = ACD_Tkr1Energy45 = 0.;

  ACD_Tkr1_VetoSigmaHit = maxSigma;
  ACD_Tkr1_VetoSigmaGap = maxSigma;

  // LOOP over AcdTkrHitPoca & get least distance sutff
  // Note that the Poca are sorted.  
  // Once we have filled all the variables we can split
  const Event::AcdTkrHitPoca* tile_vetoPoca = 0;
  const Event::AcdTkrHitPoca* ribbon_vetoPoca = 0;
  const Event::AcdTkrGapPoca* gap_vetoPoca = 0;


  float best_tileVetoHitSq(maxSigmaSq);
  float best_ribbonVetoHitSq(maxSigmaSq);
  float best_gapVetoSq(maxSigmaSq);
  float bestCornerGapMeasure(maxDist);
  
  const Event::AcdTkrAssocCol& assocs = pACD->getTkrAssocCol();
  for (Event::AcdTkrAssocCol::const_iterator itrAssoc = assocs.begin(); 
       itrAssoc != assocs.end(); itrAssoc++ ) {
    
    const Event::AcdTkrAssoc* anAssoc = (*itrAssoc);
    int trackIndex = anAssoc->getTrackIndex();
    if ( trackIndex < 0 ) continue;
    bool isBestTrack = trackIndex == 0;
    bool isUpGoing = anAssoc->getUpward();
    
    float track_tileBestVetoHitSq(maxSigmaSq);
    float track_ribbonBestVetoHitSq(maxSigmaSq);
    float track_gapBestVetoSq(maxSigmaSq);

    const Event::AcdTkrHitPoca* track_tileVetoPoca(0);
    const Event::AcdTkrHitPoca* track_ribbonVetoPoca(0);

    int nHitPoca = anAssoc->nHitPoca();    
    for ( int iHitPoca(0); iHitPoca < nHitPoca; iHitPoca++ ) {
      
      const Event::AcdTkrHitPoca* aPoca = anAssoc->getHitPoca(iHitPoca);
      if ( aPoca == 0 ) continue;
      
      // make sure there is associated signal
      if ( aPoca->mipsPmtA() < 0.001 && aPoca->mipsPmtB() < 0.001 ) continue;
            
      // get the id
      idents::AcdId theId = aPoca->getId();
      
      // We could check this earlier, but it helps to have it here, because this makes it easier to 
      // add variables associated with the down going end
      if ( isUpGoing ) {
	// Fill variables for all tracks
	float testHitSigma2 = aPoca->vetoSigma2();
	if ( theId.tile() ) {
	  if ( testHitSigma2  <  track_tileBestVetoHitSq) {
	    track_tileVetoPoca = aPoca;
	    track_tileBestVetoHitSq = testHitSigma2;
	  }	  
	} else if ( theId.ribbon() ) {
	  if ( testHitSigma2 <  track_ribbonBestVetoHitSq) {
	    track_ribbonVetoPoca = aPoca;
	    track_ribbonBestVetoHitSq = testHitSigma2;
	  }
	}      

	// Fill special vars for best track
	if ( isBestTrack ) {
	  // check to see if it is in xx deg cone around track
	  if ( theId.tile() ) {
	    if ( aPoca->getArcLength() > 0.1 ) {
	      static const float tan30 = 5.77350288616910401e-01;
	      static const float tan15 = 2.67949200239410490e-01;
	      float tanAngle = ( -1 * aPoca->getDoca() )/ aPoca->getArcLength();
	      float tileEng = hitMap[theId]->tileEnergy();
	      if ( tanAngle < 1. ) {
		ACD_Tkr1Energy45 += tileEng;
		if ( tanAngle < tan30 ) {
		  ACD_Tkr1Energy30 += tileEng;
		  if ( tanAngle < tan15 ) {
		    ACD_Tkr1Energy15 += tileEng;
		  }
		}
	      }
	    }
	  }
	}	
      }
    }

    const Event::AcdTkrGapPoca* aGap = anAssoc->getGapPoca();
    if ( aGap != 0 ) {      
      switch ( aGap->getId().gapType() )  {
      case 1: //AcdRecon::X_RibbonSide: 
      case 2: //AcdRecon::Y_RibbonSide:
      case 3: //AcdRecon::Y_RibbonTop:
      case 4: //AcdRecon::SideCornerEdge:
	track_gapBestVetoSq = aGap->vetoSigma2();
	break;
      default:
	aGap = 0;
	break;
      }
    }
      
    float cornerDoca = anAssoc->getCornerDoca();
    float cornerGapMeasure = cornerDoca > 0 ? fabs(cornerDoca) : fabs(0.2 * cornerDoca);
    if ( isUpGoing ) {
      if ( cornerGapMeasure < bestCornerGapMeasure ) {
	ACD_Corner_DOCA = cornerDoca;
	bestCornerGapMeasure = cornerGapMeasure;
      }
    }

    // Fill vars for best track
    if ( isBestTrack ) {
      
      if ( isUpGoing ) {
	ACD_Tkr1Corner_DOCA = cornerDoca;
	ACD_Tkr1_VetoSigmaHit = track_tileBestVetoHitSq < track_ribbonBestVetoHitSq ? 
	  sqrt(track_tileBestVetoHitSq) : sqrt(track_ribbonBestVetoHitSq);
      }
 
      if ( track_tileVetoPoca != 0 && isUpGoing ) {
	ACD_Tkr1ActiveDist = track_tileVetoPoca->getDoca();
	ACD_Tkr1ActiveDist_Err = track_tileVetoPoca->getDocaErrProj();
	ACD_Tkr1ActiveDist_Energy = hitMap[track_tileVetoPoca->getId()]->tileEnergy();
	ACD_Tkr1ActiveDist_Arc = track_tileVetoPoca->getArcLength();
      }
      
      if ( track_ribbonVetoPoca != 0 && isUpGoing ) {
	ACD_Tkr1_ribbon_ActiveDist =  track_ribbonVetoPoca->getDoca();
	ACD_Tkr1_ribbon_ActiveDist_Err = track_ribbonVetoPoca->getDocaErrProj();
	ACD_Tkr1_ribbon_ActiveLength = track_ribbonVetoPoca->getLocalY();
	ACD_Tkr1_ribbon_EnergyPmtA = hitMap[track_ribbonVetoPoca->getId()]->ribbonEnergy(Event::AcdHit::A);
	ACD_Tkr1_ribbon_EnergyPmtB = hitMap[track_ribbonVetoPoca->getId()]->ribbonEnergy(Event::AcdHit::B);
      }

      if ( aGap != 0 && isUpGoing ) {
	ACD_Tkr1_VetoSigmaGap = sqrt(track_gapBestVetoSq);
	ACD_Tkr1Ribbon_Dist = aGap->getDoca();
	ACD_Tkr1Ribbon_Dist_Err = aGap->getDocaErrProj();
	ACD_Tkr1RibbonLength = aGap->getLocalY();
      }
    }

    // Latch vars for all tracks
    if ( isUpGoing ) {
      if ( track_tileBestVetoHitSq < best_tileVetoHitSq  ) {
	best_tileVetoHitSq = track_tileBestVetoHitSq;
	tile_vetoPoca = track_tileVetoPoca;
      }
      if (track_ribbonBestVetoHitSq < best_ribbonVetoHitSq  ) {
	best_ribbonVetoHitSq = track_ribbonBestVetoHitSq;
	ribbon_vetoPoca = track_ribbonVetoPoca;
      }
      if ( track_gapBestVetoSq < best_gapVetoSq ) {
	best_gapVetoSq = track_gapBestVetoSq;
	gap_vetoPoca = aGap;
      }	    
    }
  }

  // Now fill in the values for the most likely Track-Tile Veto Poca
  if( tile_vetoPoca != 0 ) {
    ACD_ActiveDist3D = tile_vetoPoca->getDoca(); 
    ACD_ActiveDist3D_Err = tile_vetoPoca->getDocaErrProj();	    
    idents::AcdId theId = tile_vetoPoca->getId();
    ACD_ActiveDist_Energy = hitMap[tile_vetoPoca->getId()]->tileEnergy();
    ACD_ActiveDist_TrackNum = tile_vetoPoca->trackIndex(); // Index starts from 0
  }
  
  // Now fill in the values for the most likely Track-Ribbon Veto Poca
  if ( ribbon_vetoPoca != 0 ) {
    ACD_ribbon_ActiveDist =  ribbon_vetoPoca->getDoca();
    ACD_ribbon_ActiveDist_Err = ribbon_vetoPoca->getDocaErrProj();
    ACD_ribbon_ActiveLength = ribbon_vetoPoca->getLocalY();
    ACD_ribbon_EnergyPmtA = hitMap[ribbon_vetoPoca->getId()]->ribbonEnergy(Event::AcdHit::A);
    ACD_ribbon_EnergyPmtB = hitMap[ribbon_vetoPoca->getId()]->ribbonEnergy(Event::AcdHit::B);
  }

  if ( gap_vetoPoca != 0 ) {
    ACD_TkrRibbon_Dist = gap_vetoPoca->getDoca();
    ACD_TkrRibbon_Dist_Err = gap_vetoPoca->getDocaErrProj();
    ACD_TkrRibbonLength = gap_vetoPoca->getLocalY();    
  }

  ACD_Tkr_VetoSigmaHit = best_tileVetoHitSq < best_ribbonVetoHitSq ?
    sqrt(best_tileVetoHitSq) : sqrt(best_ribbonVetoHitSq);
  ACD_Tkr_VetoSigmaGap = sqrt(best_gapVetoSq);
    
  return sc;
}


void Acd2ValsTool::loadTopologyVars(const Event::AcdReconV2& pACD) {

  const Event::AcdEventTopology& evtTopo = pACD.getEventTopology();

  ACD_Tile_Count = evtTopo.getTileCount();
  ACD_Ribbon_Count = evtTopo.getRibbonCount();

  ACD_Veto_Count = evtTopo.getTileVeto();
  ACD_Veto_Faces = evtTopo.getNSidesVeto();

  ACD_Total_Energy = evtTopo.getTileEnergy();
  ACD_Total_Ribbon_Energy = evtTopo.getRibbonEnergy();
 
  ACD_tileTopCount = evtTopo.getVetoCountTop();
  ACD_tileCount0 = evtTopo.getVetoCountSideRow(0);
  ACD_tileCount1 = evtTopo.getVetoCountSideRow(1);
  ACD_tileCount2 = evtTopo.getVetoCountSideRow(2);
  ACD_tileCount3 = evtTopo.getVetoCountSideRow(3);
  
  ACD_countRow3Readout = evtTopo.getTileCountSideRow(3);
  
  ACD_energyTop = evtTopo.getTileEnergyTop();
  ACD_energyRow0 = evtTopo.getTileEnergySideRow(0);
  ACD_energyRow1 = evtTopo.getTileEnergySideRow(1);
  ACD_energyRow2 = evtTopo.getTileEnergySideRow(2);
  ACD_energyRow3 = evtTopo.getTileEnergySideRow(3);
}


void Acd2ValsTool::reconId(const Event::AcdReconV2 *pACD)  {

  std::vector<Event::AcdTkrHitPoca*> bestTrackUp;
  idents::AcdId resetId;
  resetId.na(1);
  ACD_TileIdRecon = resetId.id();
  ACD_RibbonIdRecon = resetId.id();
  
  // loop over the AcdTrackIntersections
  const Event::AcdTkrAssocCol& trackAssocCol = pACD->getTkrAssocCol();
  for (Event::AcdTkrAssocCol::const_iterator itrAssoc = trackAssocCol.begin();
       itrAssoc != trackAssocCol.end(); itrAssoc++ ) {
    
    // could be vertex (-1) or best track (0)            
    int trackIndex = (*itrAssoc)->getTrackIndex();
    if (trackIndex > 0 ) continue;
    Event::AcdTkrHitPoca* poca = const_cast<Event::AcdTkrHitPoca*>((*itrAssoc)->getHitPoca());
    if ( poca == 0 ) continue;	
    bool upGoing = (*itrAssoc)->getUpward();
    
    if (trackIndex != 0) {
      continue;
    }
    if (! upGoing) {
      continue;
    }
    bestTrackUp.push_back(poca);
  }
  
  setId(bestTrackUp,ACD_TileIdRecon, false);
  setId(bestTrackUp,ACD_RibbonIdRecon, true);
}

void Acd2ValsTool::setId(const std::vector<Event::AcdTkrHitPoca*> bestTrackUp,
			 unsigned int &retId, bool findRibbon) {

    idents::AcdId tUpId;
    tUpId.na(1);
    
    findId(bestTrackUp, tUpId, findRibbon);
    if ( !tUpId.na() ) {
      retId = tUpId.id();
    }
}

void Acd2ValsTool::findId(const std::vector<Event::AcdTkrHitPoca*>& vec, 
			  idents::AcdId &retId, bool findRibbon) {
  
  //Point tilePos, ribbonPos;
  int tileIndex = -1, ribIndex = -1;
  
  retId.na(1);
  if (vec.size() <= 0) {
    return;
  }
  
  // If there is only one TkrIntesection object, then we can return the one id we have
  if (vec.size() == 1) {
    if ( (!findRibbon) && (vec[0]->getId().tile()) )
      retId = vec[0]->getId();
    else if ( (findRibbon) && (vec[0]->getId().ribbon()) ) 
      retId = vec[0]->getId();
    return;
  }
  
    // otherwise, loop over the remaining TkrIntesection objects
  unsigned int ind=0;
  for ( std::vector<Event::AcdTkrHitPoca*>::const_iterator itr = vec.begin(); 
	itr != vec.end(); itr++ ) {
    idents::AcdId id = (*itr)->getId();
    if ( (!findRibbon) && (id.tile()) ) {
      if (tileIndex < 0) {  // haven't seen another tile yet
	tileIndex = ind;
	retId = id;
      } else { // there was another tile found already
	// use Z coordinates to choose if one of the found tiles is a "top" tile
	// chose the greater Z value
	if ( (retId.top()) || (id.top()) ) {
	  if (vec[tileIndex]->getGlobalPosition().z() < (*itr)->getGlobalPosition().z()) {
	    retId = (*itr)->getId();
	    tileIndex = ind;
	  }
	}   
	// assume side tiles do not overlap in Gleam, so no worries 
	// about handling that case now
	// right now we just pick up the tile we find first
      }
    } else if (findRibbon && id.ribbon() ) { // ribbon
      if (ribIndex < 0) { // first ribbon intersection found
	ribIndex = ind;
	retId = id;
      } 
      // don't worry about overlapping ribbons, 
      // just pick up the first ribbon we find
    }    
    ind++;
  }
  return;
}

