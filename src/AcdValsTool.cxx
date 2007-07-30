
/** @file AcdValsTool.cxx
@brief Calculates the Adc analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/AcdValsTool.cxx,v 1.41 2007/06/07 17:00:13 lsrea Exp $
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
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Digi/AcdDigi.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
// Point used by AcdDigi
//#include "CLHEP/Geometry/Point3D.h"  //<=== check this
#include "CLHEP/Geometry/Transform3D.h"
// Point used by TKR
//#include "geometry/Point.h"  <=== check this

#include "AcdUtil/AcdTileFuncs.h"
//#include "AcdRecon/AcdGap.h" (Don't include AcdRecon, use ints instead of enums)

#include <algorithm>
#include <numeric>
#include <map>

/** @class AcdValsTool
@brief Calculates Acd Values
*/

class AcdValsTool : public ValBase
{
public:

    AcdValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~AcdValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

    void tkrHitsCount();

    void setId(const std::vector<Event::AcdTkrIntersection*> vertexUp, 
                        const std::vector<Event::AcdTkrIntersection*> vertexDown,
                        const std::vector<Event::AcdTkrIntersection*> trackUp,
                        const std::vector<Event::AcdTkrIntersection*> trackDown,
                        unsigned int &retId, bool findRibbon=false);


    void findId(const std::vector<Event::AcdTkrIntersection*> vec, idents::AcdId &retId, bool findRibbon=false);

    void reconId(const Event::AcdRecon *pACD);

private:

    //Global ACDTuple Items
    unsigned int ACD_Tile_Count; 
    unsigned int ACD_Ribbon_Count;
    float ACD_Total_Energy;
    float ACD_Total_Ribbon_Energy;
    unsigned int ACD_TileIdRecon;
    unsigned int ACD_RibbonIdRecon;

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
    float ACD_TkrRibbonLength; 

    // Variables computed by taking best track w.r.t. hit tiles
    float ACD_Tkr1ActiveDist;
    float ACD_Tkr1ActiveDist_Err;
    float ACD_Tkr1ActiveDist_Energy;

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
    float ACD_Tkr1RibbonLength;

    // Variables computed by taking vertex w.r.t. hit tiles
    float ACD_VtxActiveDist;
    float ACD_VtxActiveDist_Energy;
    float ACD_VtxActiveDist_Down;
    float ACD_VtxActiveDist_EnergyDown;

    // Variables about number of ACD tiles by row
    unsigned int ACD_tileTopCount;
    unsigned int ACD_tileCount0;
    unsigned int ACD_tileCount1;
    unsigned int ACD_tileCount2;
    unsigned int ACD_tileCount3;

    // Services
    IGlastDetSvc *m_detSvc;

    // Algorithm parameters
    double m_vetoThresholdMeV;
};


/** @page anatup_vars
    @section adcvalstool AdCValsTool Variables
    Notes
    - Default Doca/ActiveDistance is -2000.
    - Active distance is negative if a track is outside a tile, 
    positive if inside.

<table>
<tr><th> Variable <th> Type <th> Description					

<tr><td> AcdTileCount 
<td>I<td>   Number of tiles fired
<tr><td> AcdRibbonCount 
<td>I<td>   Number of ribbons fired
<tr><td> AcdTotalEnergy 	
<td>F<td>   Total energy deposited in ACD Tiles
<tr><td> AcdRibbonEnergy 	
<td>F<td>   Total energy deposited in ACD Ribbons
<tr><td> AcdTileIdRecon
<td>I<td> Tile identifier that was pierced by the reconstructed track.  
          A value of 899 (N/A) is the default and denotes that no ACD tile was 
		  intersected by a reconstructed track.
<tr><td> AcdRibbonIdRecon (fixme)
<td>I<td> Ribbon identifier that was pierced by the reconstructed track.  
          A value of 899 (N/A) is the default and denotes that no ACD ribbon was 
		  intersected by a reconstructed track.

<tr><td> AcdActiveDist3D  (fixme)	
<td>F<td>   Largest active distance of any track to the edge of any tile 
<tr><td> AcdActiveDist3DErr (fixme)
<td>F<td>   Error on largest active distance of any track to the edge of any tile 
<tr><td> AcdActDistTileEnergy 
<td>F<td>   The deposited energy in the corresponding hit tile 

<tr><td> AcdRibbonActDist   (fixme)
<td>F<td>   Largest active distance to any ribbon 
            (considered as a straight line of no thickness) 
<tr><td> AcdRibbonActDistErr   (fixme)
<td>F<td>   Error on the smallest active distance to any ribbon 
            (considered as a straight line of no thickness) 
<tr><td> AcdRibbonActLength(fixme)
<td>F<td>   Length along ribbon where point of closest approach occured. 
            0 is center of ribbon + going towards +x or +y side of ACD
<tr><td> AcdRibbonActEnergyPmtA   (fixme)
<td>F<td>   The deposited energy in the A PMT of the corresponding hit ribbon
<tr><td> AcdRibbonActEnergyPmtB   (fixme)
<td>F<td>   The deposited energy in the B PMT of the corresponding hit ribbon

<tr><td> AcdCornerDoca 
<td>F<td>   Minimum Distance of Closest Approach of best track to the corner side gaps 
<tr><td> AcdTkrHoleDist (fixme)
<td>F<td>   Minimum Distance of Closest Approach of any track to any of the tile screw holes
<tr><td> AcdTkrRibbonDist
<td>F<td>   Minimum Distance of Closest Approach of any track to any ribbons that cover gaps
<tr><td> AcdTkrRibbonLength (fixme)
<td>F<td>   Length along ribbon where point of closest approach occured.
            0 is center of ribbon + going towards +x or +y side of ACD.

<tr><td>AcdTkr1ActiveDist 
<td>F<td>   Largest active distance from  track 1 to the edge of any tile
<tr><td>AcdTkr1ActiveDistErr
<td>F<td>   Error on largest active distance from track 1 to the edge of any tile
<tr><td>AcdTkr1ActDistTileEnergy
<td>F<td>   The deposited energy in the corresponding hit tile

<tr><td> AcdTkr1RibbonActDist   
<td>F<td>   Largest active distance to any ribbon 
            (considered as a straight line of no thickness) 
<tr><td> AcdTkr1RibbonActDistErr   
<td>F<td>   Error on the smallest active distance to any ribbon 
            (considered as a straight line of no thickness) 
<tr><td> AcdTkr1RibbonActLength
<td>F<td>   Length along ribbon where point of closest approach occured. 
            0 is center of ribbon + going towards +x or +y side of ACD
<tr><td> AcdTkr1RibbonActEnergyPmtA   
<td>F<td>   The deposited energy in the A PMT of the corresponding hit ribbon
<tr><td> AcdTkr1RibbonActEnergyPmtB   
<td>F<td>   The deposited energy in the A PMT of the corresponding hit ribbon

<tr><td> AcdTkr1CornerDoca 
<td>F<td>   Minimum Distance of Closest Approach of best track to the corner side gaps 
<tr><td> AcdTkr1HoleDist
<td>F<td>   Minimum Distance of Closest Approach to best track to any of the tile screw holes
<tr><td> AcdTkr1RibbonDist
<td>F<td>   Minimum Distance of Closest Approach to best track to any ribbons that cover gaps 
<tr><td> AcdTkr1RibbonLength (fixme)
<td>F<td>   Minimum Distance of Closest Approach to best track to any ribbons that cover gaps 

<tr><td>AcdVtxActiveDist
<td>F<td>   Largest active distance from vertex extrapolation to the edge of any tile
<tr><td>AcdVtxActiveDist_Down
<td>F<td>   Largest active distance from vertex extrapolation to the edge of any tile, down going side of tracks
<tr><td>AcdVtxActDistTileEnergy
<td>F<td>   The deposited energy in the corresponding hit tile
<tr><td>AcdVtxActDistTileEnergy_Down
<td>F<td>   The deposited energy in the corresponding hit tile, down going side of tracks

<tr><td> AcdNoSideRow[0...3] 
<td>F<td>   Hit Tile count for side row [0...3] 

</table>
    */

// predicate to identify top, (row -1) or  side  (row 0-2)
namespace {
  class acd_row { 
  public:
    acd_row(int row):m_row(row){}
    bool operator() ( std::pair<idents::AcdId ,double> entry){
      return m_row==-1? entry.first.face() == 0 : entry.first.row()==m_row;
    }
    int m_row;
  };
} 

// Static factory for instantiation of algtool objects
static ToolFactory<AcdValsTool> s_factory;
const IToolFactory& AcdValsToolFactory = s_factory;

// Standard Constructor
AcdValsTool::AcdValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 

    // in MeV
    declareProperty("VetoThresholdMeV", m_vetoThresholdMeV=0.0);
}

StatusCode AcdValsTool::initialize()
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
            StatusCode sc = m_detSvc->getNumericConstByName("acd.vetoThreshold", &m_vetoThresholdMeV);
            if (sc.isFailure()) {
                log << MSG::INFO << "Unable to retrieve threshold, setting the value to 0.4 MeV" << endreq;
                m_vetoThresholdMeV = 0.4;

            }
        }
    } else {
        return StatusCode::FAILURE;
    }

    // load up the map
    addItem("AcdTileCount",    &ACD_Tile_Count);
    addItem("AcdRibbonCount", &ACD_Ribbon_Count);
    addItem("AcdTotalEnergy", &ACD_Total_Energy);
    addItem("AcdRibbonEnergy", &ACD_Total_Ribbon_Energy);
    addItem("AcdTileIdRecon", &ACD_TileIdRecon);
    addItem("AcdRibbonIdRecon", &ACD_RibbonIdRecon);

    addItem("AcdActiveDist3D",   &ACD_ActiveDist3D);
    addItem("AcdActiveDist3DErr",   &ACD_ActiveDist3D_Err);
    addItem("AcdActDistTileEnergy",   &ACD_ActiveDist_Energy);

    addItem("AcdRibbonActDist", &ACD_ribbon_ActiveDist);
    addItem("AcdRibbonActDistErr", &ACD_ribbon_ActiveDist_Err);
    addItem("AcdRibbonActLength", &ACD_ribbon_ActiveLength);
    addItem("AcdRibbonActEnergyPmtA", &ACD_ribbon_EnergyPmtA);
    addItem("AcdRibbonActEnergyPmtB", &ACD_ribbon_EnergyPmtB);

    addItem("AcdCornerDoca",    &ACD_Corner_DOCA);
    addItem("AcdTkrHoleDist",    &ACD_TkrHole_Dist);
    addItem("AcdTkrRibbonDist",    &ACD_TkrRibbon_Dist);
    addItem("AcdTkrRibbonLength",    &ACD_TkrRibbonLength);

    addItem("AcdTkr1ActiveDist", &ACD_Tkr1ActiveDist);
    addItem("AcdTkr1ActiveDistErr", &ACD_Tkr1ActiveDist_Err);
    addItem("AcdTkr1ActDistTileEnergy", &ACD_Tkr1ActiveDist_Energy);
    
    addItem("AcdTkr1RibbonActDist", &ACD_Tkr1_ribbon_ActiveDist);
    addItem("AcdTkr1RibbonActDistErr", &ACD_Tkr1_ribbon_ActiveDist_Err);
    addItem("AcdTkr1RibbonActLength", &ACD_Tkr1_ribbon_ActiveLength);
    addItem("AcdTkr1RibbonActEnergyPmtA", &ACD_Tkr1_ribbon_EnergyPmtA);
    addItem("AcdTkr1RibbonActEnergyPmtB", &ACD_Tkr1_ribbon_EnergyPmtB);

    addItem("AcdTkr1CornerDoca",    &ACD_Tkr1Corner_DOCA);
    addItem("AcdTkr1HoleDist",    &ACD_Tkr1Hole_Dist);
    addItem("AcdTkr1RibbonDist",    &ACD_Tkr1Ribbon_Dist);
    addItem("AcdTkr1RibbonLength",    &ACD_Tkr1RibbonLength);    

    addItem("AcdVtxActiveDist", &ACD_VtxActiveDist);
    addItem("AcdVtxActDistTileEnergy", &ACD_VtxActiveDist_Energy);
    addItem("AcdVtxActiveDist_Down", &ACD_VtxActiveDist_Down);
    addItem("AcdVtxActDistTileEnergy_Down", &ACD_VtxActiveDist_EnergyDown);

    addItem("AcdNoTop",        &ACD_tileTopCount);
    addItem("AcdNoSideRow0",   &ACD_tileCount0);
    addItem("AcdNoSideRow1",   &ACD_tileCount1);
    addItem("AcdNoSideRow2",   &ACD_tileCount2);   
    addItem("AcdNoSideRow3",   &ACD_tileCount3);   

    zeroVals();

    return sc;
}


StatusCode AcdValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::AcdRecon>           pACD(m_pEventSvc,EventModel::AcdRecon::Event);

    // Recover Track associated info. (not currently used 
    //SmartDataPtr<Event::TkrFitTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    //SmartDataPtr<Event::TkrVertexCol>     pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    // Recover pointer to ACD info  

    //Make sure we have valid ACD data
    if (pACD)
    {
        reconId(pACD);

        // Make a map relating AcdId to energy in the tile
        std::map<idents::AcdId, double> tileEnergyIdMap;
	std::map<idents::AcdId, std::pair<double,double> > ribbonEnergyIdMap;

	const Event::AcdHitCol& hitCol = pACD->getAcdHitCol();
	int nHit = hitCol.size();

	static const float MeVMipTile10 = 1.9;
	static const float MeVMipTile12 = 2.28;
	static const float MeVMipRibbon = 0.5;
	
	ACD_Total_Energy = 0.;
	ACD_Total_Ribbon_Energy = 0.;
	ACD_Tile_Count = 0;
	ACD_Ribbon_Count = 0;
	
	// Loop over the hits and fill the maps
	for (int iHit(0); iHit < nHit; iHit++ ){
	  const Event::AcdHit* aHit = hitCol[iHit];
	  const idents::AcdId& id = aHit->getAcdId();
	  if ( id.na() ) {
	    continue;
	  } else if ( id.ribbon() ) {
	    float MeV_A = aHit->mips(Event::AcdHit::A) * MeVMipRibbon;
	    float MeV_B = aHit->mips(Event::AcdHit::B) * MeVMipRibbon;
	    ACD_Total_Ribbon_Energy += (MeV_A + MeV_B);
	    ribbonEnergyIdMap[id] = std::pair<double,double>(MeV_A,MeV_B);
	    ACD_Ribbon_Count++;
	  } else {
	    // check for 10mm v. 12mm tiles
	    float MeVMip = id.top() && id.row() == 2 ? MeVMipTile12 : MeVMipTile10;
	    float MeV = aHit->mips() * MeVMip;
	    ACD_Total_Energy += MeV;
	    tileEnergyIdMap[id] = MeV;
	    ACD_Tile_Count++;
	  }	  
	}
	
	// Reset variables for loop over all tracks
	ACD_ActiveDist3D = ACD_ribbon_ActiveDist = -2000.;
	ACD_ActiveDist_Energy = ACD_ribbon_EnergyPmtA = ACD_ribbon_EnergyPmtB = 0.;
	ACD_ActiveDist3D_Err = ACD_ribbon_ActiveDist_Err = -1.;	
	ACD_Corner_DOCA = ACD_TkrRibbon_Dist = ACD_TkrHole_Dist = -2000.;
	ACD_ribbon_ActiveLength = ACD_TkrRibbonLength = -10000.;

	// Reset variables for best track
	ACD_Tkr1ActiveDist = ACD_Tkr1_ribbon_ActiveDist = -2000.;
	ACD_Tkr1ActiveDist_Energy = ACD_Tkr1_ribbon_EnergyPmtA = ACD_Tkr1_ribbon_EnergyPmtB = 0.;
	ACD_Tkr1ActiveDist_Err = ACD_Tkr1_ribbon_ActiveDist_Err = -1.;
	ACD_Tkr1Corner_DOCA = ACD_Tkr1Ribbon_Dist = ACD_Tkr1Hole_Dist = -2000.;
	ACD_Tkr1_ribbon_ActiveLength = ACD_Tkr1RibbonLength = -10000.;

	// Reset vertex variables
	ACD_VtxActiveDist = ACD_VtxActiveDist_Down = -2000.;
	ACD_VtxActiveDist_Energy = ACD_VtxActiveDist_EnergyDown = 0.;

	// loop over AcdTkrHitPoca & get least distance sutff
	// Note that the Poca are sorted.  Once we have filled all the variables we can split
	const Event::AcdTkrHitPocaCol& pocas = pACD->getAcdTkrHitPocaCol();
	for ( Event::AcdTkrHitPocaCol::const_iterator itrPoca = pocas.begin(); 
	      itrPoca != pocas.end(); itrPoca++ ) {

	  const Event::AcdTkrHitPoca* aPoca = (*itrPoca);
	  // check to see if this is the vertex (-1) or the best track (0) and direction
	  bool isVertex = aPoca->trackIndex() == -1;
	  bool isBestTrack = aPoca->trackIndex() == 0;
	  bool isUpGoing = aPoca->getArcLength() > 0.;

	  bool doneHole = false;
	  double holeDoca(0.), holeDocaError(0.); 
	  int iHole = -1;

	  bool donePlaneError = false;
	  double planeError = 0.;

	  // get the id
	  idents::AcdId theId = aPoca->getId();
	  
	  if ( isUpGoing ) {
	    // Fill variables for all tracks
	    if ( theId.tile() ) {
	      if ( ACD_ActiveDist3D < -1999.99 ) {
		if ( ! doneHole ) {
		  doneHole = true;
		  AcdTileUtil::tileScrewHoleDoca(aPoca->getId(),aPoca->getActiveX(),aPoca->getActiveY(),
						 aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),aPoca->getLocalXYCov(),
						 holeDoca,holeDocaError,iHole);
		}
		if ( ! donePlaneError ) {
		  donePlaneError = true;
		  AcdTileUtil::planeErrorProjection(aPoca->getActiveX(),aPoca->getActiveY(),
						    aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),planeError);  
		}
		ACD_ActiveDist3D = aPoca->getDoca();
		ACD_ActiveDist3D_Err = planeError;
		ACD_ActiveDist_Energy = tileEnergyIdMap[theId];
		ACD_TkrHole_Dist = holeDoca;
	      } 
	    } else if ( theId.ribbon() ) {
	      if ( ACD_ribbon_ActiveDist < -1999.99 ) {
		if ( ! donePlaneError ) {
		  donePlaneError = true;
		  AcdTileUtil::planeErrorProjection(aPoca->getActiveX(),aPoca->getActiveY(),
						    aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),planeError);  
		}
		ACD_ribbon_ActiveDist =  aPoca->getDoca();
		ACD_ribbon_ActiveDist_Err = planeError;
		ACD_ribbon_ActiveLength = aPoca->getActiveY();
		ACD_ribbon_EnergyPmtA = ribbonEnergyIdMap[theId].first;
		ACD_ribbon_EnergyPmtB = ribbonEnergyIdMap[theId].second;
	      }
	    }
	    // Fill variables for best track, if appropiate
	    if ( isBestTrack ) {
	      // check tile or ribbon
	      if ( theId.tile() ) {
		// check to see if vars already filled
		if ( ACD_Tkr1ActiveDist < -1999.99 ) {
		  if ( ! doneHole ) {
		    doneHole = true;
		    AcdTileUtil::tileScrewHoleDoca(aPoca->getId(),aPoca->getActiveX(),aPoca->getActiveY(),
						   aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),aPoca->getLocalXYCov(),
						   holeDoca,holeDocaError,iHole);
		  }
		  if ( ! donePlaneError ) {
		    donePlaneError = true;
		    AcdTileUtil::planeErrorProjection(aPoca->getActiveX(),aPoca->getActiveY(),
						      aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),planeError);  
		  }
		  ACD_Tkr1ActiveDist = aPoca->getDoca();
		  ACD_Tkr1ActiveDist_Err = planeError;
		  ACD_Tkr1ActiveDist_Energy = tileEnergyIdMap[theId];
		  ACD_Tkr1Hole_Dist = holeDoca;
		} 
	      } else if ( theId.ribbon() ) {
		// check to see if vars already filled
		if ( ACD_Tkr1_ribbon_ActiveDist < -1999.99 ) {
		  if ( ! donePlaneError ) {
		    donePlaneError = true;
		    AcdTileUtil::planeErrorProjection(aPoca->getActiveX(),aPoca->getActiveY(),
						      aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),planeError);  
		  }
		  ACD_Tkr1_ribbon_ActiveDist =  aPoca->getDoca();
		  ACD_Tkr1_ribbon_ActiveDist_Err = planeError;
		  ACD_Tkr1_ribbon_ActiveLength = aPoca->getActiveY();
		  ACD_Tkr1_ribbon_EnergyPmtA = ribbonEnergyIdMap[theId].first;
		  ACD_Tkr1_ribbon_EnergyPmtB = ribbonEnergyIdMap[theId].second;
		}
	      }
	      // Of fill variables for vertex, if appropraite
	    } else if ( isVertex ) {
	      // only fill this one for tiles
	      if ( theId.tile() ) {
		// check to see if vars already filled
		if ( ACD_VtxActiveDist < -1999.99 ) {
		  ACD_VtxActiveDist = aPoca->getDoca();
		  ACD_VtxActiveDist_Energy = tileEnergyIdMap[theId];
		}
	      }
	    }
	  } else {
	    // Down going, only fill vertex variables
	    if ( isVertex ) {
	      // only fill this one for tiles
	      if ( theId.tile() ) {
		// check to see if vars already filled
		if ( ACD_VtxActiveDist_Down < -1999.99 ) {
		  ACD_VtxActiveDist_Down = aPoca->getDoca();
		  ACD_VtxActiveDist_EnergyDown = tileEnergyIdMap[theId];
		}
	      }
	    }
	  }
	}
	 
	float  bestCornerGapMeasure= 2000.;

	// loop over AcdGaps & get least distance between track extrapolation & ribbon / corner gaps
	const Event::AcdTkrGapPocaCol& gaps = pACD->getAcdTkrGapPocaCol();
	for ( Event::AcdTkrGapPocaCol::const_iterator itrGap = gaps.begin(); 
	      itrGap != gaps.end(); itrGap++ ) {

	  // check to see if this is the vertex (-1) or the best track (0) and direction
	  bool isVertex = (*itrGap)->trackIndex() == -1;
	  bool isBestTrack = (*itrGap)->trackIndex() == 0;
	  bool isUpGoing = (*itrGap)->getArcLength() > 0.;
	  bool isRibbonGap = false;
	  bool isCornerGap = false;

	  // for now we ignore vertex and down going
	  if ( isVertex or ! isUpGoing ) continue;

	  // Classify gap type
	  switch ( (*itrGap)->getId().gapType() ) {
	  case 1: //AcdRecon::X_RibbonSide: 
	  case 2: //AcdRecon::Y_RibbonSide:
	  case 3: //AcdRecon::Y_RibbonTop:
	    isRibbonGap = true;
	    break;
	  case 7: //AcdRecon::CornerRay:
	    isCornerGap = true;
	    break;
	  default:
	    break;
	  }

	  float gapDoca = (*itrGap)->getDoca();
	  if ( isRibbonGap ) {
	    // Fill variables for best track
	    if ( isBestTrack ) {
	      if ( gapDoca > ACD_Tkr1Ribbon_Dist ) {
		ACD_Tkr1Ribbon_Dist = gapDoca;
		ACD_Tkr1RibbonLength = (*itrGap)->getActiveY();
	      }
	    }
	    // Fill variables for all tracks
	    if ( gapDoca > ACD_TkrRibbon_Dist ) {
	      ACD_TkrRibbon_Dist = gapDoca;
	      ACD_TkrRibbonLength = (*itrGap)->getActiveY();
	    }
	  } else if ( isCornerGap ) {
	    //if ( (*itrGap)->getArcLength() < 0. ) continue;
	    if ( isBestTrack ) {
	      if ( fabs( gapDoca ) < fabs(ACD_Tkr1Corner_DOCA) ) {
		ACD_Tkr1Corner_DOCA = gapDoca;
	      }
	    }
	    float cornerGapMeasure = gapDoca > 0 ? gapDoca : gapDoca / -5.;
	    if ( cornerGapMeasure < bestCornerGapMeasure ) {
	      ACD_Corner_DOCA = gapDoca;
	      bestCornerGapMeasure = cornerGapMeasure;
	    }
	  }
	}    


	// use acd_row predicate to count number of tiles per side row
	unsigned int m_acd_tileCount[5] = {0,0,0,0,0};
	for( int row = -1; row<=3; row++ ){ 
	  m_acd_tileCount[row+1] = std::count_if(tileEnergyIdMap.begin(), tileEnergyIdMap.end(), acd_row(row) );
	}
	ACD_tileTopCount = m_acd_tileCount[0];
        ACD_tileCount0 = m_acd_tileCount[1];
        ACD_tileCount1 = m_acd_tileCount[2];
        ACD_tileCount2 = m_acd_tileCount[3];       
        ACD_tileCount3 = m_acd_tileCount[4];  

    } else {
        return StatusCode::FAILURE;
    }




    return sc;
}


void AcdValsTool::reconId(const Event::AcdRecon *pACD) {
    std::vector<Event::AcdTkrIntersection*> bestTrackUp, bestTrackDown, vertexUp, vertexDown;

    Event::AcdTkrIntersectionCol::const_iterator itrTrackIntersect; 

    idents::AcdId resetId;
    resetId.na(1);
    ACD_TileIdRecon = resetId.id();
    ACD_RibbonIdRecon = resetId.id();

    // loop over the AcdTrackIntersections
    const Event::AcdTkrIntersectionCol& trackIntersectCol = pACD->getAcdTkrIntersectionCol();
    for ( itrTrackIntersect = trackIntersectCol.begin(); 
        itrTrackIntersect != trackIntersectCol.end(); itrTrackIntersect++ ) {

            if ((*itrTrackIntersect)->getTrackIndex() > 0 ) continue;  // could be vertex (-1) or best track (0)
           
            double arcLen = (*itrTrackIntersect)->getArcLengthToIntersection();

            if (((*itrTrackIntersect)->getTrackIndex() == 0)){ // tracks
                if (arcLen > 0) // upward track
                    bestTrackUp.push_back(*itrTrackIntersect);
                else  // downward track
                    bestTrackDown.push_back(*itrTrackIntersect);
            } else if ( ((*itrTrackIntersect)->getTrackIndex() == -1)) { // vertices 
                if (arcLen > 0)
                    vertexUp.push_back(*itrTrackIntersect);
                else 
                    vertexDown.push_back(*itrTrackIntersect);
            }
        }


        setId(vertexUp, vertexDown, bestTrackUp, bestTrackDown, ACD_TileIdRecon, false);
        setId(vertexUp, vertexDown, bestTrackUp, bestTrackDown, ACD_RibbonIdRecon, true);

        bestTrackUp.clear();
        bestTrackDown.clear();
        vertexUp.clear();
        vertexDown.clear();
}


void AcdValsTool::findId(const std::vector<Event::AcdTkrIntersection*> vec, idents::AcdId &retId, bool findRibbon) {

    //Point tilePos, ribbonPos;
    int tileIndex = -1, ribIndex = -1;

    //idents::AcdId tileId;  //, ribId;
    //tileI d.na(1);
    //ribId.na(1);
    retId.na(1);
    if (vec.size() <= 0) {
        //retId = tileId;
        return;
    }

    //if ( (!findRibbon) && (vec[0]->getTileId().tile()) {
    //    tileIndex = 0;
    //    tileId = vec[0]->getTileId();
    //} else {
    //    ribIndex = 0;
    //    ribId = vec[0]->getTileId();
   // }

    // If there is only one TkrIntesection object, then we can return the one id we have
    if (vec.size() == 1) {
        if ( (!findRibbon) && (vec[0]->getTileId().tile()) )
            retId = vec[0]->getTileId();
        else if ( (findRibbon) && (vec[0]->getTileId().ribbon()) ) 
            retId = vec[0]->getTileId();
        return;
    }

    // otherwise, loop over the remaining TkrIntesection objects
    unsigned int ind=0;
    Event::AcdTkrIntersectionCol::const_iterator itrTrackIntersect; 
    for ( itrTrackIntersect = vec.begin(); 
        itrTrackIntersect != vec.end(); itrTrackIntersect++ ) {
            idents::AcdId id = (*itrTrackIntersect)->getTileId();
            if ( (!findRibbon) && (id.tile()) ) {
                if (tileIndex < 0) {  // haven't seen another tile yet
                    tileIndex = ind;
                    retId = id;
                } else { // there was another tile found already
                    // use Z coordinates to choose if one of the found tiles is a "top" tile
                    // chose the greater Z value
                    if ( (retId.top()) || (id.top()) ) {
                        if (vec[tileIndex]->getGlobalPosition().z() < (*itrTrackIntersect)->getGlobalPosition().z()) {
                            retId = (*itrTrackIntersect)->getTileId();
                            tileIndex = ind;
                        }
                    }   
                    // assume side tiles do not overlap in Gleam, so no worries about handling that case right now
                    // right now we just pick up the tile we find first
                }
            } else if (findRibbon && id.ribbon() ) { // ribbon
                if (ribIndex < 0) { // first ribbon intersection found
                    ribIndex = ind;
                    retId = id;
                } 
                // don't worry about overlapping ribbons, just pick up the first ribbon we find
            }

            ind++;
        }


        return;
}

void AcdValsTool::setId(const std::vector<Event::AcdTkrIntersection*> vertexUp, 
                        const std::vector<Event::AcdTkrIntersection*> vertexDown,
                        const std::vector<Event::AcdTkrIntersection*> bestTrackUp,
                        const std::vector<Event::AcdTkrIntersection*> bestTrackDown,
                        unsigned int &retId, bool findRibbon) {
                            
    idents::AcdId vUpId, vDownId, tUpId, tDownId;
    vUpId.na(1);
    vDownId.na(1);
    tUpId.na(1);
    tDownId.na(1);

    findId(vertexUp, vUpId, findRibbon);
    findId(vertexDown, vDownId, findRibbon);
    findId(bestTrackUp, tUpId, findRibbon);
    findId(bestTrackDown, tDownId, findRibbon);

     // Prefer vertices over tracks
     // up versus down
     bool found = false;
     if (!vUpId.na()) {
         retId = vUpId.id();
         found = true;
      } else if (!tUpId.na()) { 
         retId = tUpId.id();
         found = true;
      }

      // no upward intersections found - check downward
      if (!found) {
          if (!vDownId.na()) 
              retId = vDownId.id();
          else if (!tDownId.na()) 
               retId = tDownId.id();
      }
}
