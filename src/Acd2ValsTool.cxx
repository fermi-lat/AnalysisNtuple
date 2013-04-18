
/** @file AcdValsTool.cxx
@brief Calculates the Adc analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/Acd2ValsTool.cxx,v 1.12 2013/04/18 20:21:19 echarles Exp $
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
#include "Event/Recon/AcdRecon/AcdAssoc.h"


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
  int          ACD_Trigger_Veto; 
  unsigned int ACD_Ribbon_Count;
  unsigned int ACD_Veto_Count; 
  unsigned int ACD_Veto_Faces;
  float        ACD_Total_Tile_Energy;
  float        ACD_Total_Ribbon_Energy;
  // ADW: More global items
  float        ACD_Tile_Energy;
  float        ACD_Ribbon_Energy;
  float        ACD_Ghost_Tile_Energy;
  float        ACD_Ghost_Ribbon_Energy;
  float        ACD_Trigger_Tile_Energy;
  float        ACD_Trigger_Ribbon_Energy;

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
  // ADW: Trigger energy in cone
  float ACD_Tkr1TriggerEnergy15;
  float ACD_Tkr1TriggerEnergy30;
  float ACD_Tkr1TriggerEnergy45;
  
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
    
  // New stuff EAC & ADW
  float ACD_Tkr1_VetoSigmaHit;
  float ACD_Tkr1_VetoSigmaGap;
  float ACD_Tkr1_VetoSigmaMip;
  float ACD_Tkr1_VetoSigmaProp;
  float ACD_Tkr1_VetoSigmaProj;
  unsigned int  ACD_Tkr1_TriggerVeto;
  float ACD_Tkr_VetoSigmaHit;
  float ACD_Tkr_VetoSigmaGap;
  float ACD_Tkr_VetoSigmaMip;
  float ACD_Tkr_VetoSigmaProp;
  float ACD_Tkr_VetoSigmaProj;
  unsigned int  ACD_Tkr_TriggerVeto;

  // Variables computed for best CAL cluster
  float ACD_Cal1_VetoSigmaHit;
  float ACD_Cal1_VetoSigmaMip;
  float ACD_Cal1_VetoSigmaProp;
  float ACD_Cal1_VetoSigmaProj;

  // Variables computed by taking best cluster w.r.t. hit tiles
  float ACD_Cal1Energy15;
  float ACD_Cal1Energy30;
  float ACD_Cal1Energy45;
  float ACD_Cal1TriggerEnergy15;
  float ACD_Cal1TriggerEnergy30;
  float ACD_Cal1TriggerEnergy45;
  
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
  
  // variables for CR active distance
  float ACD_CR_ActiveDist3D;
  float ACD_CR_ActiveDist_Energy;
  int ACD_CR_ActiveDist_TrackNum;
  float ACD_CR_ribbon_ActiveDist;
  float ACD_CR_ribbon_EnergyPmtA;
  float ACD_CR_ribbon_EnergyPmtB;
  float ACD_CR1_ActiveDist;
  float ACD_CR1_ribbon_ActiveDist;
  float ACD_CR1_ribbon_EnergyPmtA;
  float ACD_CR1_ribbon_EnergyPmtB;
  float ACD_CR1_ActiveDist_Energy;
  int ACD_CR1_ActiveDist_TrackNum;

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
- Default action is to de-ghost energies at the global level, 
but not at the individual track association level.

<table>
<tr><th> Variable <th> Type <th> Description                    
<tr><td> Acd2TriggerVeto 
<td>U<td>   Set if any ACD trigger veto fired
<tr><td> Acd2TileCount 
<td>U<td>   Number of tiles fired
<tr><td> Acd2TriggerVeto 
<td>U<td>   Trigger veto for the fast signal
<tr><td> Acd2RibbonCount 
<td>U<td>   Number of ribbons fired
<tr><td> Acd2VetoCount 
<td>U<td>   Total number of vetoes fired in ACD
<tr><td> Acd2VetoFaces
<td>U<td>   Number of ACD faces with tile vetos
<tr><td> Acd2NumTileVetoTop
<td>U<td>   Number of top tiles with a veto
<tr><td> Acd2NumTileVetoSideRow[0...3] 
<td>U<td>   Number of tiles in side row [0...3] with a veto
<tr><td> Acd2NumTotalTileHitRow3
<td>U<td>   Number of tiles with a hit for side row 3, no veto threshold applied (includes ghosts)

<tr><td> Acd2TotalTileEnergy    
<td>F<td>   Total energy deposited in ACD Tiles (includes ghosts)
<tr><td> Acd2TotalRibbonEnergy   
<td>F<td>   Total energy deposited in ACD Ribbons (includes ghosts)
<tr><td> Acd2TotalTileEnergyTop
<td>F<td>   Total energy deposited in top tiles (includes ghosts)
<tr><td> Acd2TotalTileEnergyRow[0...3]
<td>F<td>   Total energy deposited in the tiles in side row 0 to 3 (includes ghosts)

<tr><td> Acd2TileEnergy    
<td>F<td>   Energy (de-ghosted) deposited in ACD Tiles
<tr><td> Acd2RibbonEnergy   
<td>F<td>   Energy (de-ghosted) deposited in ACD Ribbons 
<tr><td> Acd2GhostTileEnergy    
<td>F<td>   Energy deposited in ACD Tiles from hits tagged as ghosts
<tr><td> Acd2GhostRibbonEnergy   
<td>F<td>   Energy deposited in ACD Ribbons from hits tagged as ghosts
<tr><td> Acd2TriggerTileEnergy    
<td>F<td>   Energy deposited in ACD Tiles with trigger veto asserted
<tr><td> Acd2TriggerRibbonEnergy   
<td>F<td>   Energy deposited in ACD Ribbons with trigger veto asserted

<tr><td> Acd2TileIdRecon
<td>U<td> Tile identifier that was pierced by the reconstructed track.  
A value of 899 (N/A) is the default and denotes that no ACD tile was 
intersected by a reconstructed track.
<tr><td> Acd2RibbonIdRecon
<td>U<td> Ribbon identifier that was pierced by the reconstructed track.  
A value of 899 (N/A) is the default and denotes that no ACD ribbon was 
intersected by a reconstructed track.


<tr><td> Acd2TileActDist3D   
<td>F<td>   Tile Active Distance most likely to give a veto.  
Corresponds to Act. Dist. that is greater than a set 
an energy dep. min. distance and has the largest pulse height
<tr><td> Acd2TileActDist3DErr
<td>F<td>   Error on most likely veto active distance of any track 
to the edge of any tile 
<tr><td> Acd2TileActDistEnergy 
<td>F<td>   The deposited energy (not de-ghosted) in the corresponding hit tile 
<tr><td> Acd2TileActDistTrackNum
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
<td>F<td>   The deposited energy (not de-ghosted) in the A PMT of the corresponding hit ribbon
<tr><td> Acd2RibbonActEnergyPmtB
<td>F<td>   The deposited energy (not de-ghosted) in the B PMT of the corresponding hit ribbon
<tr><td> Acd2CornerDoca 
<td>F<td>   Minimum Distance of Closest Approach of a track to the corner side gaps
This variable is signed to match the direction of the overlaps in the tiles. 
The gap appears larger for tracks coming from the + side than the - side.  
<tr><td> Acd2TkrHoleDist
<td>F<td>   FIXME: Minimum Distance of Closest Approach of any track to any of the tile screw holes
<tr><td> Acd2TkrRibbonDist
<td>F<td>   Minimum Distance of Closest Approach of any track to any ribbons that cover gaps
<tr><td> Acd2TkrRibbonDistErr
<td>F<td>   Error On Minimum Distance of Closest Approach of any track to any ribbons that cover gaps
<tr><td> Acd2TkrRibbonLength
<td>F<td>   Length along ribbon where point of closest approach occured.
0 is center of ribbon + going towards +x or +y side of ACD.


<tr><td> Acd2Tkr1TileActDist 
<td>F<td>   Largest active distance from  track 1 to the edge of any tile
<tr><td> Acd2Tkr1TileActDistErr
<td>F<td>   Error on largest active distance from track 1 to the edge of any tile
<tr><td> Acd2Tkr1TileActDistEnergy
<td>F<td>   The deposited energy (not de-ghosted) in the corresponding hit tile
<tr><td> Acd2Tkr1TileActDistArc
<td>F<td>   Length from head of track to point where active distances was calculated

<tr><td> Acd2Tkr1RibbonActDist   
<td>F<td>   Largest active distance to any ribbon 
(considered as a straight line at the center of the ribbon) 
<tr><td> Acd2Tkr1RibbonActDistErr   
<td>F<td>   Error on the smallest active distance to any ribbon 
<tr><td> Acd2Tkr1RibbonActLength
<td>F<td>   Length along ribbon where point of closest approach occured. 
0 is center of ribbon + going towards +x or +y side of ACD
<tr><td> Acd2Tkr1RibbonActDistEnergyPmtA   
<td>F<td>   The deposited energy (not de-ghosted) in the A PMT of the corresponding hit ribbon
<tr><td> Acd2Tkr1RibbonActDistEnergyPmtB   
<td>F<td>   The deposited energy (not de-ghosted) in the A PMT of the corresponding hit ribbon

<tr><td> Acd2Tkr1Energy15
<td>F<td>   Energy (de-ghosted) in 15 deg. cone ahead of Track
<tr><td> Acd2Tkr1Energy30
<td>F<td>   Energy (de-ghosted) in 30 deg. cone ahead of Track
<tr><td> Acd2Tkr1Energy45
<td>F<td>   Energy (de-ghosted) in 45 deg. cone ahead of Track
<tr><td> Acd2Tkr1TriggerEnergy15
<td>F<td>   Trigger energy in 15 deg. cone ahead of Track
<tr><td> Acd2Tkr1TriggerEnergy30
<td>F<td>   Trigger energy in 30 deg. cone ahead of Track
<tr><td> Acd2Tkr1TriggerEnergy45
<td>F<td>   Trigger energy in 45 deg. cone ahead of Track

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

<tr><td> Acd2Tkr1VetoSigmaMip
<td>F<td>Number of sigmas less than an expected mip for signal in tile or ribbon most likely to veto the best track.
<tr><td> Acd2Tkr1VetoSigmaProp
<td>F<td>Number of sigmas track propagation is away from tile or ribbon most likely to veto the best track.
<tr><td> Acd2Tkr1VetoSigmaProj
<td>F<td>Number of sigmas track projection is away from tile or ribbon most likely to veto the best track.
<tr><td> Acd2Tkr1VetoSigmaHit
<td>F<td>Number of sigmas less than an expected mip for signal combined with the
number of sigmas track propagation is away from tile or ribbon most likely to veto the best track.
<tr><td> Acd2Tkr1VetoSigmaGap
<td>F<td>Number of sigmas track propagation is away from clostest GAP in ACD for the best track.
<tr><td> Acd2Tkr1TriggerVeto
<td>U<td>Veto trigger for the tile or ribbon most likely to veto the best track.

<tr><td> Acd2TkrVetoSigmaMip
<td>F<td>Number of sigmas less than an expected mip for signal in tile or ribbon most likely to veto any track.
<tr><td> Acd2TkrVetoSigmaProp
<td>F<td>Number of sigmas track propagation is away from tile or ribbon most likely to veto any track.
<tr><td> Acd2TkrVetoSigmaProj
<td>F<td>Number of sigmas track projection is away from tile or ribbon most likely to veto any track.
<tr><td> Acd2TkrVetoSigmaHit
<td>F<td>Number of sigmas less than an expected mip for signal combined with the
number of sigmas track propagation is away from tile or ribbon most likely to veto any track.
<tr><td> Acd2TkrVetoSigmaGap
<td>F<td>Number of sigmas track propagation is away from closest GAP in ACD for any track.
<tr><td> Acd2Tkr1TriggerVeto
<td>U<td>Veto trigger for the tile or ribbon most likely to veto any track.

<tr><td> Acd2Cal1VetoSigmaMip
<td>F<td>Number of sigmas less than an expected mip for signal in tile or ribbon most likely to veto the first 
CAL cluster.
<tr><td> Acd2Cal1VetoSigmaProp
<td>F<td>Number of sigmas track propagation is away from tile or ribbon most likely to veto the first CAL cluster.
<tr><td> Acd2Cal1VetoSigmaProj
<td>F<td>Number of sigmas track projection is away from tile or ribbon most likely to veto the first CAL cluster.
<tr><td> Acd2Cal1VetoSigmaHit
<td>F<td>Number of sigmas less than an expected mip for signal combined with the
number of sigmas track propagation is away from tile or ribbon most likely to veto the first CAL cluster.

<tr><td> Acd2Cal1Energy15
<td>F<td>   Energy (de-ghosted) in 15 deg. cone ahead of CAL cluster
<tr><td> Acd2Cal1Energy30
<td>F<td>   Energy (de-ghosted) in 30 deg. cone ahead of CAL cluster
<tr><td> Acd2Cal1Energy45
<td>F<td>   Energy (de-ghosted) in 45 deg. cone ahead of CAL cluster
<tr><td> Acd2Cal1TriggerEnergy15
<td>F<td>   Trigger energy in 15 deg. cone ahead of CAL cluster
<tr><td> Acd2Cal1TriggerEnergy30
<td>F<td>   Trigger energy in 30 deg. cone ahead of CAL cluster
<tr><td> Acd2Cal1TriggerEnergy45
<td>F<td>   Trigger energy in 45 deg. cone ahead of CAL cluster

</table>
*/


// Static factory for instantiation of algtool objects
//static ToolFactory<Acd2ValsTool> s_factory;
//const IToolFactory& Acd2ValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(Acd2ValsTool);

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
  addItem(m_prefix + "TriggerVeto",    &ACD_Trigger_Veto);
  addItem(m_prefix + "TileCount",    &ACD_Tile_Count);
  addItem(m_prefix + "TriggerVeto",    &ACD_Trigger_Veto);
  addItem(m_prefix + "RibbonCount", &ACD_Ribbon_Count);
  addItem(m_prefix + "VetoCount",    &ACD_Veto_Count);
  addItem(m_prefix + "VetoFaces", &ACD_Veto_Faces);

  addItem(m_prefix + "NumTileVetoTop",        &ACD_tileTopCount);
  addItem(m_prefix + "NumTileVetoSideRow0",   &ACD_tileCount0);
  addItem(m_prefix + "NumTileVetoSideRow1",   &ACD_tileCount1);
  addItem(m_prefix + "NumTileVetoSideRow2",   &ACD_tileCount2);   
  addItem(m_prefix + "NumTileVetoSideRow3",   &ACD_tileCount3);
  addItem(m_prefix + "NumTotalTileHitRow3",   &ACD_countRow3Readout);

  addItem(m_prefix + "TotalTileEnergy", &ACD_Total_Tile_Energy);
  addItem(m_prefix + "TotalRibbonEnergy", &ACD_Total_Ribbon_Energy);
  addItem(m_prefix + "TotalTileEnergyTop",      & ACD_energyTop);
  addItem(m_prefix + "TotalTileEnergyRow0",     & ACD_energyRow0);
  addItem(m_prefix + "TotalTileEnergyRow1",     & ACD_energyRow1);
  addItem(m_prefix + "TotalTileEnergyRow2",     & ACD_energyRow2);
  addItem(m_prefix + "TotalTileEnergyRow3",     & ACD_energyRow3);


  addItem(m_prefix + "TileEnergy", &ACD_Tile_Energy);
  addItem(m_prefix + "RibbonEnergy", &ACD_Ribbon_Energy);
  addItem(m_prefix + "GhostTileEnergy", &ACD_Ghost_Tile_Energy);
  addItem(m_prefix + "GhostRibbonEnergy", &ACD_Ghost_Ribbon_Energy);
  addItem(m_prefix + "TriggerTileEnergy", &ACD_Trigger_Tile_Energy, true);
  addItem(m_prefix + "TriggerRibbonEnergy", &ACD_Trigger_Ribbon_Energy);
  addItem(m_prefix + "TileIdRecon", &ACD_TileIdRecon);
  addItem(m_prefix + "RibbonIdRecon", &ACD_RibbonIdRecon);
  
  addItem(m_prefix + "TileActDist3D",       &ACD_ActiveDist3D);
  addItem(m_prefix + "TileActDist3DErr",    &ACD_ActiveDist3D_Err);
  addItem(m_prefix + "TileActDistEnergy",   &ACD_ActiveDist_Energy);
  addItem(m_prefix + "TileActDistTrackNum", &ACD_ActiveDist_TrackNum);
  
  addItem(m_prefix + "RibbonActDist",       &ACD_ribbon_ActiveDist);
  addItem(m_prefix + "RibbonActDistErr",    &ACD_ribbon_ActiveDist_Err);
  addItem(m_prefix + "RibbonActLength",     &ACD_ribbon_ActiveLength);
  addItem(m_prefix + "RibbonActEnergyPmtA", &ACD_ribbon_EnergyPmtA);
  addItem(m_prefix + "RibbonActEnergyPmtB", &ACD_ribbon_EnergyPmtB);

  addItem(m_prefix + "CornerDoca",         &ACD_Corner_DOCA);
  addItem(m_prefix + "TkrHoleDist",        &ACD_TkrHole_Dist);
  addItem(m_prefix + "TkrRibbonDist",      &ACD_TkrRibbon_Dist);
  addItem(m_prefix + "TkrRibbonDistErr",   &ACD_TkrRibbon_Dist_Err);
  addItem(m_prefix + "TkrRibbonLength",    &ACD_TkrRibbonLength);
  
  addItem(m_prefix + "Tkr1TileActDist",       &ACD_Tkr1ActiveDist);
  addItem(m_prefix + "Tkr1TileActDistErr",    &ACD_Tkr1ActiveDist_Err);
  addItem(m_prefix + "Tkr1TileActDistArc",    &ACD_Tkr1ActiveDist_Arc);
  addItem(m_prefix + "Tkr1TileActDistEnergy", &ACD_Tkr1ActiveDist_Energy);

  addItem(m_prefix + "Tkr1RibbonActDist",       &ACD_Tkr1_ribbon_ActiveDist);
  addItem(m_prefix + "Tkr1RibbonActDistErr",    &ACD_Tkr1_ribbon_ActiveDist_Err);
  addItem(m_prefix + "Tkr1RibbonActLength",     &ACD_Tkr1_ribbon_ActiveLength);
  addItem(m_prefix + "Tkr1RibbonActEnergyPmtA", &ACD_Tkr1_ribbon_EnergyPmtA);
  addItem(m_prefix + "Tkr1RibbonActEnergyPmtB", &ACD_Tkr1_ribbon_EnergyPmtB);
  
  addItem(m_prefix + "Tkr1Energy15",        &ACD_Tkr1Energy15);
  addItem(m_prefix + "Tkr1Energy30",        &ACD_Tkr1Energy30, true);
  addItem(m_prefix + "Tkr1Energy45",        &ACD_Tkr1Energy45);
  addItem(m_prefix + "Tkr1TriggerEnergy15", &ACD_Tkr1TriggerEnergy15);
  addItem(m_prefix + "Tkr1TriggerEnergy30", &ACD_Tkr1TriggerEnergy30);
  addItem(m_prefix + "Tkr1TriggerEnergy45", &ACD_Tkr1TriggerEnergy45);
  
  addItem(m_prefix + "Tkr1CornerDoca",      &ACD_Tkr1Corner_DOCA, true);
  addItem(m_prefix + "Tkr1HoleDist",        &ACD_Tkr1Hole_Dist);
  addItem(m_prefix + "Tkr1RibbonDist",      &ACD_Tkr1Ribbon_Dist);
  addItem(m_prefix + "Tkr1RibbonDistErr",   &ACD_Tkr1Ribbon_Dist_Err);
  addItem(m_prefix + "Tkr1RibbonLength",    &ACD_Tkr1RibbonLength);    
  
  addItem(m_prefix + "Tkr1VetoSigmaHit",  & ACD_Tkr1_VetoSigmaHit);
  addItem(m_prefix + "Tkr1VetoSigmaGap",  & ACD_Tkr1_VetoSigmaGap,  true);
  addItem(m_prefix + "Tkr1VetoSigmaMip",  & ACD_Tkr1_VetoSigmaMip,  true);
  addItem(m_prefix + "Tkr1VetoSigmaProp", & ACD_Tkr1_VetoSigmaProp);
  addItem(m_prefix + "Tkr1VetoSigmaProj", & ACD_Tkr1_VetoSigmaProj, true);
  addItem(m_prefix + "Tkr1TriggerVeto",   & ACD_Tkr1_TriggerVeto);
  addItem(m_prefix + "TkrVetoSigmaHit",   & ACD_Tkr_VetoSigmaHit);
  addItem(m_prefix + "TkrVetoSigmaGap",   & ACD_Tkr_VetoSigmaGap);   
  addItem(m_prefix + "TkrVetoSigmaMip",   & ACD_Tkr_VetoSigmaMip);
  addItem(m_prefix + "TkrVetoSigmaProp",  & ACD_Tkr_VetoSigmaProp);
  addItem(m_prefix + "TkrVetoSigmaProj",  & ACD_Tkr_VetoSigmaProj);
  addItem(m_prefix + "TkrTriggerVeto",    & ACD_Tkr_TriggerVeto);

  addItem(m_prefix + "Cal1VetoSigmaHit", & ACD_Cal1_VetoSigmaHit);
  addItem(m_prefix + "Cal1VetoSigmaMip", & ACD_Cal1_VetoSigmaMip);
  addItem(m_prefix + "Cal1VetoSigmaProp", & ACD_Cal1_VetoSigmaProp);
  addItem(m_prefix + "Cal1VetoSigmaProj", & ACD_Cal1_VetoSigmaProj);

  addItem(m_prefix + "Cal1Energy15", &ACD_Cal1Energy15, true);
  addItem(m_prefix + "Cal1Energy30", &ACD_Cal1Energy30);
  addItem(m_prefix + "Cal1Energy45", &ACD_Cal1Energy45);
  addItem(m_prefix + "Cal1TriggerEnergy15", &ACD_Cal1TriggerEnergy15);
  addItem(m_prefix + "Cal1TriggerEnergy30", &ACD_Cal1TriggerEnergy30);
  addItem(m_prefix + "Cal1TriggerEnergy45", &ACD_Cal1TriggerEnergy45);
 
  addItem(m_prefix + "CRActiveDist3D",   &ACD_CR_ActiveDist3D);
  addItem(m_prefix + "CRActDistTileEnergy",   &ACD_CR_ActiveDist_Energy);
  addItem(m_prefix + "CRActDistTrackNum", &ACD_CR_ActiveDist_TrackNum);
  addItem(m_prefix + "CRRibbonActiveDist", &ACD_CR_ribbon_ActiveDist);
  addItem(m_prefix + "CRRibbonActEnergyPmtA", &ACD_CR_ribbon_EnergyPmtA);
  addItem(m_prefix + "CRRibbonActEnergyPmtB", &ACD_CR_ribbon_EnergyPmtB);
  addItem(m_prefix + "CR1ActiveDist",   &ACD_CR1_ActiveDist);
  addItem(m_prefix + "CR1ActDistTileEnergy",   &ACD_CR1_ActiveDist_Energy);
  addItem(m_prefix + "CR1ActDistTrackNum", &ACD_CR1_ActiveDist_TrackNum);
  addItem(m_prefix + "CR1RibbonActiveDist", &ACD_CR1_ribbon_ActiveDist);
  addItem(m_prefix + "CR1RibbonActEnergyPmtA", &ACD_CR1_ribbon_EnergyPmtA);
  addItem(m_prefix + "CR1RibbonActEnergyPmtB", &ACD_CR1_ribbon_EnergyPmtB);

 
  zeroVals();
  ACD_ActiveDist_TrackNum = -1;
  ACD_CR_ActiveDist_TrackNum = -1; //RJ
  ACD_CR1_ActiveDist_TrackNum = -1; //RJ

  return sc;
}


StatusCode Acd2ValsTool::calculate()
{

  StatusCode sc = StatusCode::SUCCESS;
  
  // Recover pointers to ACD Recon results
  SmartDataPtr<Event::AcdReconV2> pACD(m_pEventSvc,EventModel::AcdReconV2::Event);

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

  /// ADW: Incorporated into AcdEventTopology
  // Set to zero for each event (MPR)
  //int TriggerVeto = 0;

  // Loop over the hits and fill the maps

  for (int iHit(0); iHit < nHit; iHit++ ){
    const Event::AcdHit* aHit = hitCol[iHit];
    const idents::AcdId& id = aHit->getAcdId();
  
    //Trigger veto for the fast signal (MPR).
    //int triggerVetoBit_0 = ((aHit->getFlags(Event::AcdHit::A) >> 1) & 0x1);
    //int triggerVetoBit_1 = ((aHit->getFlags(Event::AcdHit::B) >> 1) & 0x1);
    //if (triggerVetoBit_0 == 1 || triggerVetoBit_1 == 1){
    //  TriggerVeto = 1;
    //}
    // 
    //ACD_Trigger_Veto = TriggerVeto;
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
  
  ACD_CR_ActiveDist3D = ACD_CR_ribbon_ActiveDist = -2000.;  //RJ
  ACD_CR1_ActiveDist = ACD_CR1_ribbon_ActiveDist = -2000.;  //RJ
  ACD_CR_ActiveDist_TrackNum = ACD_CR1_ActiveDist_TrackNum = -999; //RJ
  ACD_CR_ActiveDist_Energy = ACD_CR_ribbon_EnergyPmtA = ACD_CR_ribbon_EnergyPmtB = 0.;  //RJ
  ACD_CR1_ActiveDist_Energy = ACD_CR1_ribbon_EnergyPmtA = ACD_CR1_ribbon_EnergyPmtB = 0.;  //RJ

  // Reset variables for best track
  ACD_Tkr1ActiveDist = ACD_Tkr1_ribbon_ActiveDist = negMaxDist;
  ACD_Tkr1ActiveDist_Energy = ACD_Tkr1_ribbon_EnergyPmtA 
    = ACD_Tkr1_ribbon_EnergyPmtB = 0.;
  ACD_Tkr1ActiveDist_Err = ACD_Tkr1_ribbon_ActiveDist_Err = ACD_Tkr1Ribbon_Dist_Err = -1.;
  ACD_Tkr1ActiveDist_Arc = -1.;
  ACD_Tkr1Corner_DOCA = ACD_Tkr1Ribbon_Dist = ACD_Tkr1Hole_Dist = negMaxDist;
  ACD_Tkr1_ribbon_ActiveLength = ACD_Tkr1RibbonLength = -10000.;
  ACD_Tkr1Energy15 = ACD_Tkr1Energy30 = ACD_Tkr1Energy45 = 0.;
  ACD_Tkr1TriggerEnergy15 = ACD_Tkr1TriggerEnergy30 = ACD_Tkr1TriggerEnergy45 = 0.;

  ACD_Tkr1_VetoSigmaHit = maxSigma;
  ACD_Tkr1_VetoSigmaGap = maxSigma;
  ACD_Tkr1_VetoSigmaMip = maxSigma;
  ACD_Tkr1_VetoSigmaProp = maxSigma;
  ACD_Tkr1_VetoSigmaProj = maxSigma;

  ACD_Tkr1_TriggerVeto = 0;

  ACD_Tkr_VetoSigmaHit = maxSigma;
  ACD_Tkr_VetoSigmaGap = maxSigma;
  ACD_Tkr_VetoSigmaMip = maxSigma;
  ACD_Tkr_VetoSigmaProp = maxSigma;
  ACD_Tkr_VetoSigmaProj = maxSigma;

  ACD_Tkr_TriggerVeto = 0;

  // RJ: Find the best CR track
  std::vector<Event::TkrTrack*> trackVec = m_pTrackVec->getTrackVec();

  int iCnt = 0;
  int bestCRtkr = -9999;
  int mxHts= 0;
  for (int iCnt(0); iCnt < trackVec.size(); iCnt++ ) {
    Event::TkrTrack* myTrack = trackVec[iCnt];
    if (myTrack->getStatusBits() & Event::TkrTrack::COSMICRAY) {
      if (myTrack->getNumHits() > mxHts) {
	mxHts= myTrack->getNumHits();
	bestCRtkr= iCnt;
      }
    }
  }


  // LOOP over AcdTkrHitPoca & get least distance stuff
  // Note that the Poca are sorted.  
  // Once we have filled all the variables we can split
  const Event::AcdTkrHitPoca* tile_vetoPoca = 0;
  const Event::AcdTkrHitPoca* ribbon_vetoPoca = 0;
  const Event::AcdTkrGapPoca* gap_vetoPoca = 0;

  const Event::AcdTkrHitPoca* tile_CR_vetoPoca = 0;
  const Event::AcdTkrHitPoca* ribbon_CR_vetoPoca = 0;


  float best_tileVetoSq(maxSigmaSq);
  float best_ribbonVetoSq(maxSigmaSq);
  float best_gapVetoSq(maxSigmaSq);

  float best_CR_tileVetoSq(maxSigmaSq);
  float best_CR_ribbonVetoSq(maxSigmaSq);

  float bestCornerGapMeasure(maxDist);
  
  const Event::AcdTkrAssocCol& assocs = pACD->getTkrAssocCol();
  for (Event::AcdTkrAssocCol::const_iterator itrAssoc = assocs.begin(); 
       itrAssoc != assocs.end(); itrAssoc++ ) {
    
    const Event::AcdAssoc* anAssoc = (*itrAssoc);
    int trackIndex = anAssoc->getTrackIndex();
    if ( trackIndex < 0 ) continue;
    bool isBestTrack = trackIndex == 0;
    bool isCR = trackIndex >= 100;
    if (isCR) {
      trackIndex -= 100;
    }
    bool isBestCR = trackIndex == bestCRtkr;
    ACD_CR1_ActiveDist_TrackNum = trackIndex;
    bool isUpGoing = anAssoc->getUpward();    

    // These are the best veto values if the element is
    // a tile (not the best tile veto)
    float track_tileBestVetoSq(maxSigmaSq);
    float track_tileBestVetoHit(maxSigma);
    float track_tileBestVetoProp(maxSigma);
    float track_tileBestVetoProj(maxSigma);
    unsigned int track_tileBestTriggerVeto(0);
    // These are the best veto values if the element is
    // a ribbon (not the best ribbon veto)
    float track_ribbonBestVetoSq(maxSigmaSq);
    float track_ribbonBestVetoHit(maxSigma);
    float track_ribbonBestVetoProp(maxSigma);
    float track_ribbonBestVetoProj(maxSigma);
    float track_gapBestVetoSq(maxSigmaSq);
    unsigned int track_ribbonBestTriggerVeto(0);

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
          if ( testHitSigma2  <  track_tileBestVetoSq) {
            track_tileVetoPoca = aPoca;
            track_tileBestVetoSq = testHitSigma2;
            track_tileBestVetoHit = aPoca->vetoSigmaHit();
            track_tileBestVetoProp = aPoca->vetoSigmaProp();
            track_tileBestVetoProj = aPoca->vetoSigmaProj();
            track_tileBestTriggerVeto = hitMap[theId]->getTriggerVeto();
          }   
        } else if ( theId.ribbon() ) {
          if ( testHitSigma2 <  track_ribbonBestVetoSq) {
            track_ribbonVetoPoca = aPoca;
            track_ribbonBestVetoSq = testHitSigma2;
            track_ribbonBestVetoHit = aPoca->vetoSigmaHit();
            track_ribbonBestVetoProp = aPoca->vetoSigmaProp();
            track_ribbonBestVetoProj = aPoca->vetoSigmaProj();
            track_ribbonBestTriggerVeto = hitMap[theId]->getTriggerVeto();
          }
        }      

        // Fill special vars for best track
        if ( isBestTrack ) {
          // check to see if it is in xx deg cone around track
          // ADW: Use de-ghosted tile energy and fill trigger tile energy
          if ( theId.tile() ) {
            if ( aPoca->getArcLength() > 0.1 ) {
              static const float tan30 = 5.77350288616910401e-01;
              static const float tan15 = 2.67949200239410490e-01;
              float tanAngle = ( -1 * aPoca->getDoca() )/ aPoca->getArcLength();
              float tileEng = hitMap[theId]->tileEnergy() * ( !hitMap[theId]->getGhost() );
              float triggerTileEng = tileEng * ( hitMap[theId]->getTriggerVeto() );
              if ( tanAngle < 1. ) {
                ACD_Tkr1Energy45 += tileEng;
                ACD_Tkr1TriggerEnergy45 += triggerTileEng;
                if ( tanAngle < tan30 ) {
                  ACD_Tkr1Energy30 += tileEng;
                  ACD_Tkr1TriggerEnergy30 += triggerTileEng;
                  if ( tanAngle < tan15 ) {
                    ACD_Tkr1Energy15 += tileEng;
                    ACD_Tkr1TriggerEnergy15 += triggerTileEng;
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
       
      if (track_tileBestVetoSq < track_ribbonBestVetoSq) { 
          ACD_Tkr1_VetoSigmaHit = sqrt(track_tileBestVetoSq);
          ACD_Tkr1_VetoSigmaMip =  track_tileBestVetoHit > 0. ?
            track_tileBestVetoHit : 0.;
          ACD_Tkr1_VetoSigmaProp = track_tileBestVetoProp > 0. ?
            track_tileBestVetoProp : 0.;
          ACD_Tkr1_VetoSigmaProj = track_tileBestVetoProj > 0. ?
            track_tileBestVetoProj : 0.;
          ACD_Tkr1_TriggerVeto = track_tileBestTriggerVeto;
        } else {
          ACD_Tkr1_VetoSigmaHit = sqrt(track_ribbonBestVetoSq);
          ACD_Tkr1_VetoSigmaMip =  track_ribbonBestVetoHit > 0. ?
            track_ribbonBestVetoHit : 0.;
          ACD_Tkr1_VetoSigmaProp = track_ribbonBestVetoProp > 0. ?
            track_ribbonBestVetoProp : 0.;
          ACD_Tkr1_VetoSigmaProj = track_ribbonBestVetoProj > 0. ?
            track_ribbonBestVetoProj : 0.;
          ACD_Tkr1_TriggerVeto = track_ribbonBestTriggerVeto;
        }
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
       
    if ( isBestCR ) {
      if ( track_tileVetoPoca != 0 && isUpGoing ) {
        ACD_CR1_ActiveDist = track_tileVetoPoca->getDoca();
        ACD_CR1_ActiveDist_Energy = hitMap[track_tileVetoPoca->getId()]->tileEnergy();
      }
        
      if ( track_ribbonVetoPoca != 0 && isUpGoing ) {
        ACD_CR1_ribbon_ActiveDist =  track_ribbonVetoPoca->getDoca();
        ACD_CR1_ribbon_EnergyPmtA = hitMap[track_ribbonVetoPoca->getId()]->ribbonEnergy(Event::AcdHit::A);
        ACD_CR1_ribbon_EnergyPmtB = hitMap[track_ribbonVetoPoca->getId()]->ribbonEnergy(Event::AcdHit::B);
      }
    }

    // Latch vars for all tracks
    if ( isUpGoing ) {
      if ( isCR ) {
	if ( track_tileBestVetoSq < best_CR_tileVetoSq  ) {
	  best_CR_tileVetoSq = track_tileBestVetoSq;
	  tile_CR_vetoPoca = track_tileVetoPoca;
	}
	if (track_ribbonBestVetoSq < best_CR_ribbonVetoSq  ) {
	  best_CR_ribbonVetoSq = track_ribbonBestVetoSq;
	  ribbon_CR_vetoPoca = track_ribbonVetoPoca;
	}
      } else {
	if ( track_tileBestVetoSq < best_tileVetoSq  ) {
	  best_tileVetoSq = track_tileBestVetoSq;
	  tile_vetoPoca = track_tileVetoPoca;
	}
	if (track_ribbonBestVetoSq < best_ribbonVetoSq  ) {
	  best_ribbonVetoSq = track_ribbonBestVetoSq;
	  ribbon_vetoPoca = track_ribbonVetoPoca;
	}
	if ( track_gapBestVetoSq < best_gapVetoSq ) {
	  best_gapVetoSq = track_gapBestVetoSq;
	  gap_vetoPoca = aGap;
	}     
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


  // Now fill in the values for the most likely CR Track-Tile Veto Poca
  if( tile_CR_vetoPoca != 0 ) {
    ACD_CR_ActiveDist3D = tile_CR_vetoPoca->getDoca(); 
    idents::AcdId theId = tile_CR_vetoPoca->getId();
    ACD_CR_ActiveDist_Energy = hitMap[tile_CR_vetoPoca->getId()]->tileEnergy();
    ACD_CR_ActiveDist_TrackNum = tile_CR_vetoPoca->trackIndex(); // Index starts from 0
  }
  
  // Now fill in the values for the most likely CR Track-Ribbon Veto Poca
  if ( ribbon_CR_vetoPoca != 0 ) {
    ACD_CR_ribbon_ActiveDist =  ribbon_CR_vetoPoca->getDoca();
    ACD_CR_ribbon_EnergyPmtA = hitMap[ribbon_CR_vetoPoca->getId()]->ribbonEnergy(Event::AcdHit::A);
    ACD_CR_ribbon_EnergyPmtB = hitMap[ribbon_CR_vetoPoca->getId()]->ribbonEnergy(Event::AcdHit::B);
  }

  if ( gap_vetoPoca != 0 ) {
    ACD_TkrRibbon_Dist = gap_vetoPoca->getDoca();
    ACD_TkrRibbon_Dist_Err = gap_vetoPoca->getDocaErrProj();
    ACD_TkrRibbonLength = gap_vetoPoca->getLocalY();    
  }

  if (best_tileVetoSq < best_ribbonVetoSq) {
    ACD_Tkr_VetoSigmaHit = sqrt(best_tileVetoSq);
    if ( tile_vetoPoca != 0 ) {
      ACD_Tkr_VetoSigmaMip  =  tile_vetoPoca->vetoSigmaHit();      
      ACD_Tkr_VetoSigmaProj =  tile_vetoPoca->vetoSigmaProp();     
      ACD_Tkr_VetoSigmaProp =  tile_vetoPoca->vetoSigmaProj();     
      ACD_Tkr_TriggerVeto = hitMap[tile_vetoPoca->getId()]->getTriggerVeto();
    }
  } else {
    ACD_Tkr_VetoSigmaHit = sqrt(best_ribbonVetoSq);
    if ( ribbon_vetoPoca != 0 ) {
      ACD_Tkr_VetoSigmaMip  =  ribbon_vetoPoca->vetoSigmaHit();      
      ACD_Tkr_VetoSigmaProj =  ribbon_vetoPoca->vetoSigmaProp();     
      ACD_Tkr_VetoSigmaProp =  ribbon_vetoPoca->vetoSigmaProj();     
      ACD_Tkr_TriggerVeto = hitMap[ribbon_vetoPoca->getId()]->getTriggerVeto();
    }
  }
  //ACD_Tkr_VetoSigmaHit = best_tileVetoSq < best_ribbonVetoSq ?
  //  sqrt(best_tileVetoSq) : sqrt(best_ribbonVetoSq);

  ACD_Tkr_VetoSigmaGap = sqrt(best_gapVetoSq);

  /// ADW: 2/2/12
  /// Calorimeter Associations

  // Reset variables for first CAL cluster
  ACD_Cal1Energy15 = ACD_Cal1Energy30 = ACD_Cal1Energy45 = 0.;

  ACD_Cal1_VetoSigmaHit = maxSigma;
  ACD_Cal1_VetoSigmaMip = maxSigma;
  ACD_Cal1_VetoSigmaProp = maxSigma;
  ACD_Cal1_VetoSigmaProj = maxSigma;

  // LOOP over AcdTkrHitPoca & get least distance stuff
  // Note that the Poca are sorted.  
  // Once we have filled all the variables we can split
  tile_vetoPoca = 0;
  tile_CR_vetoPoca = 0;
  ribbon_vetoPoca = 0;
  ribbon_CR_vetoPoca = 0;

  best_tileVetoSq = maxSigmaSq;
  best_ribbonVetoSq = maxSigmaSq;
  best_CR_tileVetoSq = maxSigmaSq;
  best_CR_ribbonVetoSq = maxSigmaSq;

  const Event::AcdCalAssocCol& calAssocs = pACD->getCalAssocCol();
  for (Event::AcdCalAssocCol::const_iterator itrAssoc = calAssocs.begin(); 
       itrAssoc != calAssocs.end(); itrAssoc++ ) {

    const Event::AcdAssoc* anAssoc = (*itrAssoc);
    int clusterIndex = anAssoc->getTrackIndex();
    if ( clusterIndex < 0 ) continue;
    bool isBestCluster = clusterIndex == 0;
    bool isUpGoing = anAssoc->getUpward();
    
    // These are the best veto values if the element is
    // a tile (not the best tile veto)
    float cluster_tileBestVetoSq(maxSigmaSq);
    float cluster_tileBestVetoHit(maxSigma);
    float cluster_tileBestVetoProp(maxSigma);
    float cluster_tileBestVetoProj(maxSigma);
    // These are the best veto values if the element is
    // a ribbon (not the best ribbon veto)
    float cluster_ribbonBestVetoSq(maxSigmaSq);
    float cluster_ribbonBestVetoHit(maxSigma);
    float cluster_ribbonBestVetoProp(maxSigma);
    float cluster_ribbonBestVetoProj(maxSigma);

    const Event::AcdTkrHitPoca* cluster_tileVetoPoca(0);
    const Event::AcdTkrHitPoca* cluster_ribbonVetoPoca(0);

    // Remove this to calculate variables for all clusters
    if ( !isBestCluster ) continue;

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
          if ( testHitSigma2  <  cluster_tileBestVetoSq) {
            cluster_tileVetoPoca = aPoca;
            cluster_tileBestVetoSq = testHitSigma2;
            cluster_tileBestVetoHit = aPoca->vetoSigmaHit();
            cluster_tileBestVetoProp = aPoca->vetoSigmaProp();
            cluster_tileBestVetoProj = aPoca->vetoSigmaProj();
          }   
        } else if ( theId.ribbon() ) {
          if ( testHitSigma2 <  cluster_ribbonBestVetoSq) {
            cluster_ribbonVetoPoca = aPoca;
            cluster_ribbonBestVetoSq = testHitSigma2;
            cluster_ribbonBestVetoHit = aPoca->vetoSigmaHit();
            cluster_ribbonBestVetoProp = aPoca->vetoSigmaProp();
            cluster_ribbonBestVetoProj = aPoca->vetoSigmaProj();
          }
        }      

        // Fill special vars for best track
        if ( isBestCluster ) {
          // check to see if it is in xx deg cone around track
          if ( theId.tile() ) {
            if ( aPoca->getArcLength() > 0.1 ) {
              static const float tan30 = 5.77350288616910401e-01;
              static const float tan15 = 2.67949200239410490e-01;
              float tanAngle = ( -1 * aPoca->getDoca() )/ aPoca->getArcLength();
              float tileEng = hitMap[theId]->tileEnergy();
              float triggerTileEng = tileEng * ( hitMap[theId]->getTriggerVeto() );
              if ( tanAngle < 1. ) {
                ACD_Cal1Energy45 += tileEng;
                ACD_Cal1TriggerEnergy45 += triggerTileEng;
                if ( tanAngle < tan30 ) {
                  ACD_Cal1Energy30 += tileEng;
                  ACD_Cal1TriggerEnergy30 += triggerTileEng;
                  if ( tanAngle < tan15 ) {
                    ACD_Cal1Energy15 += tileEng;
                    ACD_Cal1TriggerEnergy15 += triggerTileEng;
                  }
                }
              }
            }
          }
        }

      }

    }

    // Fill vars for best track
    if ( isBestCluster ) {
      if ( isUpGoing ) {
        if (cluster_tileBestVetoSq < cluster_ribbonBestVetoSq) { 
          ACD_Cal1_VetoSigmaHit = sqrt(cluster_tileBestVetoSq);
          ACD_Cal1_VetoSigmaMip =  cluster_tileBestVetoHit > 0. ?
            cluster_tileBestVetoHit : 0.;
          ACD_Cal1_VetoSigmaProp = cluster_tileBestVetoProp > 0. ?
            cluster_tileBestVetoProp : 0.;
          ACD_Cal1_VetoSigmaProj = cluster_tileBestVetoProj > 0. ?
            cluster_tileBestVetoProj : 0.;
        } else {
          ACD_Cal1_VetoSigmaHit = sqrt(cluster_ribbonBestVetoSq);
          ACD_Cal1_VetoSigmaMip =  cluster_ribbonBestVetoHit > 0. ?
            cluster_ribbonBestVetoHit : 0.;
          ACD_Cal1_VetoSigmaProp = cluster_ribbonBestVetoProp > 0. ?
            cluster_ribbonBestVetoProp : 0.;
          ACD_Cal1_VetoSigmaProj = cluster_ribbonBestVetoProj > 0. ?
            cluster_ribbonBestVetoProj : 0.;
        }
      }
    }

  }    

  return sc;
}


void Acd2ValsTool::loadTopologyVars(const Event::AcdReconV2& pACD) {

  const Event::AcdEventTopology& evtTopo = pACD.getEventTopology();

  ACD_Tile_Count = evtTopo.getTileCount();
  ACD_Ribbon_Count = evtTopo.getRibbonCount();

  ACD_Veto_Count = evtTopo.getVetoCount();
  ACD_Trigger_Veto = (ACD_Veto_Count > 0);

  ACD_Veto_Faces = evtTopo.getNSidesVeto();

  ACD_Total_Tile_Energy = evtTopo.getTotalTileEnergy();
  ACD_Total_Ribbon_Energy = evtTopo.getTotalRibbonEnergy();

  ACD_Tile_Energy = evtTopo.getTileEnergy();
  ACD_Ribbon_Energy = evtTopo.getRibbonEnergy();

  ACD_Ghost_Tile_Energy = evtTopo.getGhostTileEnergy();
  ACD_Ghost_Ribbon_Energy = evtTopo.getGhostRibbonEnergy();

  ACD_Trigger_Tile_Energy = evtTopo.getTriggerTileEnergy();
  ACD_Trigger_Ribbon_Energy = evtTopo.getTriggerRibbonEnergy();
 
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

