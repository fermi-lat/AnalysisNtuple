/** @file CalValsTool.cxx
@brief Calculates the Cal analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/CalValsTool.cxx,v 1.96 2010/01/25 17:59:26 usher Exp $
*/
//#define PRE_CALMOD 1

// Include files


#include "ValBase.h"

#include "UBinterpolate.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalParams.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

#include "GaudiKernel/IToolSvc.h"
//#include "GlastSvc/Reco/IPropagatorTool.h"
//#include "GlastSvc/Reco/IPropagator.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "idents/TowerId.h" 
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/Rotation.h"

#include "TMath.h"

#include "Doca.h"

/*! @class CalValsTool
@brief calculates Cal values

@authors Bill Atwood, Leon Rochester
*/

namespace {
    double truncFraction = 0.9;
    const int _badInt = -1;
    const float _badFloat = -2.0;
}

class CalValsTool :   public ValBase
{
public:

    CalValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~CalValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

    void zeroVals();

private:

    // some pointers to services  
    /// GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc; 
    /// TkrGeometrySvc used for access to tracker geometry info
    ITkrGeometrySvc* m_tkrGeom;

    IPropagator*     m_G4PropTool;    

    IValsTool*       m_pTkrTool;

    /// some Geometry
    double m_towerPitch;
    int    m_xNum;
    int    m_yNum;
    int    m_nLayers;
    int    m_nCsI;

    /// gets the CAL info from detModel
    StatusCode getCalInfo();
    double activeDist(Point pos, int &view) const;

    /// CAL vars
    double m_calXWidth;
    double m_calYWidth;
    double m_deltaSide;
    double m_calZTop;
    double m_calZBot;

    float CAL_EnergyRaw;
    float CAL_EnergyCorr;
    float CAL_Leak_Corr;
    float CAL_Edge_Corr;
    float CAL_Total_Corr; 

    float CAL_CsI_RLn;
    float CAL_Tot_RLn;
    float CAL_Cntr_RLn; 
    float CAL_LAT_RLn; 
    float CAL_DeadTot_Rat;
    float CAL_DeadCnt_Rat; 

    float CAL_t_Pred; 
    float CAL_deltaT;
    float CAL_xEcntr;
    float CAL_yEcntr;
    float CAL_zEcntr;
    float CAL_xdir;
    float CAL_ydir;
    float CAL_zdir;
    float CAL_x0;
    float CAL_y0;
	float CAL_X0_Lyr[8];
	float CAL_Y0_Lyr[8];
    float CAL_Top_Gap_Dist;

    float CAL_xEcntr2;
    float CAL_yEcntr2;
    float CAL_zEcntr2;
    float CAL_xdir2;
    float CAL_ydir2;
    float CAL_zdir2;
    float CAL_posdir_chisq;
    float CAL_posdir_nlayers;
    float CAL_nsaturated;

    float CAL_Gap_Fraction;  
    float CAL_TwrEdgeCntr;
    float CAL_TwrEdgeTop;
    float CAL_LATEdge; 
    float CAL_EdgeEnergy;

    float CAL_Lyr0_Ratio;
    float CAL_Lyr7_Ratio;
    float CAL_BkHalf_Ratio;

    float CAL_Xtal_Ratio;
    float CAL_Xtal_maxEne; 
    float CAL_eLayer[8];
	float CAL_LightAsym[8];
    float CAL_Num_Xtals;
    float CAL_Num_Xtals_Trunc;
    float CAL_Max_Num_Xtals_In_Layer;
    float CAL_Long_Rms;
    float CAL_Trans_Rms;
    float CAL_LRms_Asym;
    float CAL_Long_Skew;
    float CAL_Long_Skew_Norm;

    float CAL_MIP_Diff; 
    float CAL_MIP_Ratio;

    // Stuff to study clustering
    float CAL_energy_uber;
    float CAL_xEcntr_uber;
    float CAL_yEcntr_uber;
    float CAL_zEcntr_uber;
    float CAL_xdir_uber;
    float CAL_ydir_uber;
    float CAL_zdir_uber;

    float CAL_xEcntr2_uber;
    float CAL_yEcntr2_uber;
    float CAL_zEcntr2_uber;
    float CAL_xdir2_uber;
    float CAL_ydir2_uber;
    float CAL_zdir2_uber;

    float CAL_num_clusters;
    float CAL_rest_energy;
    float CAL_rest_numXtals;

    // New variables for new energy correction tools
    // Full profile fit
    float CAL_cfp_energy;      // Energy from Full Profile tool
    float CAL_cfp_totChiSq;    // Total ChiSquare of fit divided by 11
    float CAL_cfp_energyUB;    // Energy from Full Profile tool unbiased
    float CAL_cfp_calEffRLn;   // Effective radiation lengths in the Cal
    float CAL_cfp_tkrRLn;   // Effective radiation lengths in the tkr
    float CAL_cfp_alpha;       // fit parameter alpha
    float CAL_cfp_tmax;        // fit parameter tmax
    float CAL_cfp_fiterrflg;   // fit error flag
    // Same variables but with fit performed with cal direction
    float CAL_cfp_calfit_energy;      // Energy from Full Profile tool
    float CAL_cfp_calfit_totChiSq;    // Total ChiSquare of fit divided by 11
    float CAL_cfp_calfit_calEffRLn;   // Effective radiation lengths in the Cal
    float CAL_cfp_calfit_alpha;       // fit parameter alpha
    float CAL_cfp_calfit_tmax;        // fit parameter tmax
    float CAL_cfp_calfit_fiterrflg;   // fit error flag

    //float CAL_lll_energy;      // Energy from the Last Layer Likelihood tool
    //float CAL_lll_energyErr;   // Chisquare from the Last Layer Likelihood
    //float CAL_tkl_energy;      // Energy from the tracker Likelihood tool
    //float CAL_tkl_energyErr;   // Chisquare from the tracker likelihood
    float CAL_LkHd_energy;     // Energy from the Likelihood tool
    float CAL_LkHd_energyErr;  // Chisquare from the likelihood
    float CAL_LkHd_energyUB;   // Energy from the Likelihood tool unbiased
    float CAL_RmsE;            // Rms of layer energies, corrected for pathlength
    //    and dropping layers with low E, normalized
    //    to the energy in these layers.
    float CAL_numLayersRms;    // Number of layers participating in CAL_RmsE

    //Calimeter items with Recon - Tracks
    float CAL_Track_DOCA;
    float CAL_Track_Angle;
    float CAL_Track_Sep;
    float CAL_track_rms;
    float CAL_track_E_rms;
    float CAL_track_rms_trunc;
    float CAL_track_E_rms_trunc;

    float CAL_rmsLayerE;
    float CAL_rmsLayerEBack;
    int   CAL_nLayersRmsBack;
    float CAL_eAveBack;
    float CAL_layer0Ratio;
    float CAL_xPosRmsLastLayer;
    float CAL_yPosRmsLastLayer;

  // Used to determine transverse size using only transverse position measurement (Philippe Bruel)
  int TSnlog;
  int TSiused[1536];
  int TSiorder[1536];
  double TSdistTL[1536];
  double TSdistT[1536];
  double TSdist[1536];
  double TSTS[1536];
  double TSenergy[1536];
  double TSefrac[1536];
  Point TSaxisP;
  Vector TSaxisV;

  float CAL_TS_CAL_TL_68;
  float CAL_TS_CAL_TL_90;
  float CAL_TS_CAL_TL_95;
  float CAL_TS_CAL_TL_99;
  float CAL_TS_CAL_TL_100;
  float CAL_TS_CAL_T_68;
  float CAL_TS_CAL_T_90;
  float CAL_TS_CAL_T_95;
  float CAL_TS_CAL_T_99;
  float CAL_TS_CAL_T_100;

  float CAL_TS_CAL2_TL_68;
  float CAL_TS_CAL2_TL_90;
  float CAL_TS_CAL2_TL_95;
  float CAL_TS_CAL2_TL_99;
  float CAL_TS_CAL2_TL_100;
  float CAL_TS_CAL2_T_68;
  float CAL_TS_CAL2_T_90;
  float CAL_TS_CAL2_T_95;
  float CAL_TS_CAL2_T_99;
  float CAL_TS_CAL2_T_100;

  float CAL_TS_TKR_TL_68;
  float CAL_TS_TKR_TL_90;
  float CAL_TS_TKR_TL_95;
  float CAL_TS_TKR_TL_99;
  float CAL_TS_TKR_TL_100;
  float CAL_TS_TKR_T_68;
  float CAL_TS_TKR_T_90;
  float CAL_TS_TKR_T_95;
  float CAL_TS_TKR_T_99;
  float CAL_TS_TKR_T_100;

  int TSfillTSdist(Event::CalXtalRecCol *pxtalrecs);
  int TSfillTS(int optts);
  double TSgetinterpolationTS(double efrac);

    UBinterpolate* m_ubInterpolate;
};

namespace {
    // this is the test distance for the CAL_EdgeEnergy variable
    double _deltaEdge = 50.0;
}

// Static factory for instantiation of algtool objects
static ToolFactory<CalValsTool> s_factory;
const IToolFactory& CalValsToolFactory = s_factory;

// Standard Constructor
CalValsTool::CalValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

/** @page anatup_vars 
    @section calvalstool CalValsTool Variables
<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> CalEnergyRaw 
<td>F<td>   Sum of the raw energies in all the crystals.  
Includes estimate of missed energy due to zero-supression.  
This replaces the variable CalEnergySum.
<tr><td> CalEnergyCorr 
<td>F<td>   Cal Energy corrected layer-by-layer for edges and leakage.  
This replaces the variable CalEneSumCorr.
<tr><td> CalLeakCorr 
<td>F<td>   Leakage correction: this is the contained fraction of the total energy 
after edge corrections.  
<tr><td> CalEdgeCorr 
<td>F<td>   Effective layer-by-layer edge correction mainly due to the gaps 
between Cal modules; multiplicative 
<tr><td> CalTotalCorr 
<td>F<td>   Global total correction. Includes effect due to dead material; 
multiplicative 
<tr><td> CalCsIRLn 
<td>F<td>   Total radiation lengths in crystals, integrated along the 
event axis (line connecting the first hit in the tracker to the CAL energy centroid) 
<tr><td> CalTotRLn 
<td>F<td>   Total radiation lengths in the CAL, integrated along the event axis. 
<tr><td> CalCntRLn 
<td>F<td>   Radiation lengths integrated along the event axis, up to energy centroid 
<tr><td> CalLATRLn 
<td>F<td>   Total radiation lengths integrated along the event axis 
(including the tracker). 
<tr><td> CalDeadTotRat 
<td>F<td>   Ratio of radiation lengths in dead material to CalTotRLn 
<tr><td> CalDeadCntRat 
<td>F<td>   Ratio of radiation lengths in dead material up to energy centroid, 
to CalCntRat 
<tr><td> CalTPred 
<td>F<td>   Model-predicted energy centroid in radiation lengths 
<tr><td> CalDeltaT 
<td>F<td>   Difference between measured and predicted energy centroids 
<tr><td> CalTwrEdge 
<td>F<td>   Distance of the entry point of the best track from the tower boundary, 
measured at the top of the CAL. 
<tr><td> CalLATEdge 
<td>F<td>   Closest distance of track 1, projected to the top of the CAL, 
to the edge of the CAL layer, taking non-square shape into account. 
This is essentially the old merit skirt variable. 
<tr><td> CalEdgeEnergy 
<td>F<td>   The sum of the raw energies in each crystal for which the energy centroid 
is within _deltaEdge (currently 50 mm) of the outside edge of one of 
the outside CAL modules.
This is an attempt at a "anti-coincidence counter" for the CAL.
<tr><td> CalTwrEdgeCntr 
<td>F<td>   Distance of the energy centroid from the nearest tower boundary.  
<tr><td> CalGapFraction 
<td>F<td>   Approximate fraction of the shower volumn which falls in inter-tower gaps. 
<tr><td> CalTrackSep 
<td>F<td>   Distance between impact points of two best tracks at CAL front face; 
zero if only one track 
<tr><td> CalTrackDoca 
<td>F<td>   Distance between the projected vertex (or track if only one track) 
and the energy centroid, evaluated at the z of the centroid. 
<tr><td> CalTrackAngle 
<td>F<td>   Angle between "gamma" direction in the tracker and direction of the CAL "track" 
<tr><td> CalELayerN, N=0,7 
<td>F<td>   Energy deposited in layer N of the CAL 
<tr><td> CalLyr0Ratio 
<td>F<td>   Ratio of CalELayer0 to CalEnergyRaw 
<tr><td> CalLyr7Ratio 
<td>F<td>   Ratio of CalELayer7 to CalEnergyRaw 
<tr><td> CalBkHalfRatio 
<td>F<td>   Ratio of total energy in back half of CAL (layers 4-7) to 
CalEnergyRaw 
<tr><td> CalNumXtals 
<td>F<td>   Number of CAL Xtals above threshold
<tr><td> CalMaxNumXtalsInLayer
<td>F<td>   Number of xtals with energy above 5 MeV in the layer
            with the highest number of xtals above 5 MeV. The cut
            is chosen to match the one used in GCRSelect
<tr><td> CalXtalsTrunc 
<td>F<td>   Number of CAL Xtals with > %1 of CalEnergyRaw 
            in first Cluster (so far, there's only one!) 
<tr><td> CalXtalRatio 
<td>F<td>   Ratio of number of Xtals with energy > 1% of CalEnergyRaw to 
total number of struck Xtals in the first Cluster  
<tr><td> CalXtalMaxEne 
<td>F<td>   Maximum energy found in a single Xtal
<tr><td> CalLongRms 
<td>F<td>   rms of the average of the 1st and 3rd shower moments. 
Indicates the length of the measured shower along the shower axis. 
<tr><td> CalLRmsAsym 
<td>F<td>   The asymetry of the 1st and 3rd shower moments.  
This should be close to zero. Because of ordering of moments it is slightly ... (??)
<tr><td> CalTransRms 
<td>F<td>   rms of transverse position measurements.
<tr><td> CalMIPDiff 
<td>F<td>   Difference between measured energy and that expected 
from a minimum-ionizing particle 
<tr><td> CalMIPRatio 
<td>F<td>   Ratio of measured energy to that expected from a 
minimum-ionizing particle 
<tr><td> Cal[X/Y/Z]Ecntr 
<td>F<td>   Energy centroid in [x/y/z]
<tr><td> Cal[X/Y/Z]Dir 
<td>F<td>   [x/y/z] direction cosine of CAL "track" 
<tr><td> Cal[X/Y]0 
<td>F<td>   [x/y] position of where the first track hits TOP of Calorimeter
<tr><td> CalTopGapDist 
<td>F<td>   Distance of TOP position from an intermodule gap
<tr><td> CalTrkXtalRms
<td>F<td>   For this and the next three variables, a measure of the rms spread
of the crystals in the calorimeter around the projection of the first
track. They differ by the weighting and crystal count. This one is 
the rms of the DOCAs of all crystals with repect to the projected track
<tr><td> CalTrkXtalRmsE
<td>F<td>   Same as the previous, but weighted by the deposited energy in 
each crystal
<tr><td> CalTrkXtalRmsTrunc
<td>F<td>   Same as the first, but excludes the crystals with the largest DOCAs.
Default is to remove 10% of the crystals, with the number of crystals
included rounded to the nearest integer
<tr><td> CalTrkXtalRmsETrunc
<td>F<td>   Same as the previous, but weighted by the deposited energy in
each crystal
<tr><td> CalRmsLayerE
<td>F<td>   Rms of deposited energy (normalized to rad. lengths traversed, 
and to average energy deposited) for all
the layers with > 5% of deposited raw energy and at least 0.5 rad lengths
of predicted CsI traversed. 
<tr><td> CalRmsLayerEBack
<td>F<td>   Same as above, with layer 0 left out of the calculation
<tr><td> CalNLayersRmsBacj
<td>I<td>    Number of layers used in the calculation of CalRmsLayerEBack above
<tr><td> CalEAveBack
<td>F<td>   Average of normalized eDep excluding layer 0, cuts as for CalRmsLayerE above
<tr><td> CalLayer0Ratio
<td>F<td>   Ratio of layer0 normalized eDep to the average of the remaining layers, 
cuts as for CalRmsLayerE above
<tr><td> Cal[X/Y]PosRmsLL
<td>F<td>   Energy Weighted Rms of the hit crystals in the last layer
(layer 7). X is measured across the width of the crystals,
Y across the length.
<tr><td> CalCfpEnergy
<td>F<td> energy reported by the full profile fit (using the tracker direction if there is one, the cal direction otherwise)
<tr><td> CalCfpChiSq
<td>F<td> chi squared reported by the full profile fit (using the tracker direction if there is one, the cal direction otherwise)
<tr><td> CalCfpEffRLn
<td>F<td> effective amount of X0 reported by the full profile fit (using the tracker direction if there is one, the cal direction otherwise)
<tr><td> CalCfpTkrRLn
<td>F<td> amount of X0 in the tracker reported by the full profile fit (using the tracker direction if there is one, the cal direction otherwise)
<tr><td> CalCfpAlpha
<td>F<td> alpha (first parameter) reported by the full profile fit (using the tracker direction if there is one, the cal direction otherwise)
<tr><td> CalCfpTmax
<td>F<td> tmax (second parameter = position of maximum of shower in X0) reported by the full profile fit (using the tracker direction if there is one, the cal direction otherwise)
<tr><td> CalCfpFitErrFlg
<td>F<td> fit flag reported by the full profile fit (using the tracker direction if there is one, the cal direction otherwise)
<tr><td> CalCfpCalEnergy
<td>F<td> energy reported by the full profile fit (using the cal direction)
<tr><td> CalCfpCalChiSq
<td>F<td> chi squared reported by the full profile fit (using the cal direction)
<tr><td> CalCfpCalEffRLn
<td>F<td> effective amount of X0 reported by the full profile fit (using the cal direction)
<tr><td> CalCfpCalAlpha
<td>F<td> alpha (first parameter) reported by the full profile fit (using the cal direction)
<tr><td> CalCfpCalTmax
<td>F<td> tmax (second parameter = position of maximum of shower in X0) reported by the full profile fit (using the cal direction)
<tr><td> CalCfpCalFitErrFlg
<td>F<td> fit flag reported by the full profile fit (using the cal direction)
<tr><td> Cal[X/Y/Z]Ecntr2 
<td>F<td>   Energy centroid in [x/y/z] determined using only the crystal transverse information
<tr><td> Cal[X/Y/Z]Dir2 
<td>F<td>   [x/y/z] direction cosine of CAL "track" determined using only the crystal transverse information
<tr><td> CalPosDirChisq
<td>F<td> chisquared of fit determination of position and direction using only the crystal transverse information
<tr><td> CalPosDirNLayers
<td>F<td> number of layers used during fit determination of position and direction using only the crystal transverse information
<tr><td> CalNSaturated
<td>F<td> number of saturated crystals
<tr><td> CalTrSizeCalT[68,90,95,99,100] 
<td>F<td> Transverse size of shower as function of energy fraction : Cal =  using cal position and direction; T = using only crystal transverse information
<tr><td> CalTrSizeCalTL[68,90,95,99,100] 
<td>F<td> Transverse size of shower as function of energy fraction : Cal =  using cal position and direction; TL = using crystal transverse and longitudinal information
<tr><td> CalTrSizeCal2T[68,90,95,99,100] 
<td>F<td> Transverse size of shower as function of energy fraction : Cal2 =  using cal2 position and direction; T = using only crystal transverse information
<tr><td> CalTrSizeTkrT[68,90,95,99,100] 
<td>F<td> Transverse size of shower as function of energy fraction : Tkr =  using tkr position and direction; T = using only crystal transverse information
<tr><td> CalTrSizeTkrTL[68,90,95,99,100] 
<td>F<td> Transverse size of shower as function of energy fraction : Tkr =  using tkr position and direction; TL = using crystal transverse and longitudinal information
</table>

*/


StatusCode CalValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the services

    if( serviceLocator() ) {

        // find GlastDevSvc service
        if (service("GlastDetSvc", m_detSvc, true).isFailure()){
            log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
            return StatusCode::FAILURE;
        }

        m_detSvc->getNumericConstByName("CALnLayer", &m_nLayers);
        m_detSvc->getNumericConstByName("nCsIPerLayer", &m_nCsI);

        // find TkrGeometrySvc service
        if (service("TkrGeometrySvc", m_tkrGeom, true).isFailure()){
            log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
            return StatusCode::FAILURE;
        }

        IToolSvc* toolSvc = 0;
        if(service("ToolSvc", toolSvc, true).isFailure()) {
            log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
            return StatusCode::FAILURE;
        }

        if(!toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)) {
            log << MSG::ERROR << "Couldn't find the G4PropagationTool!" << endreq;
            return StatusCode::FAILURE;
        }


        m_towerPitch = m_tkrGeom->towerPitch();
        m_xNum = m_tkrGeom->numXTowers();
        m_yNum = m_tkrGeom->numYTowers();

        IToolSvc* pToolSvc = 0; 
        sc = service("ToolSvc", pToolSvc, true);
        if ( !sc.isSuccess() ) {
            log << MSG::ERROR << "Can't find ToolSvc, will quit now" << endreq;
            return StatusCode::FAILURE;
        }

        m_pTkrTool = 0;
        sc = pToolSvc->retrieveTool("TkrValsTool", m_pTkrTool);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "Unable to find tool: " "TkrValsTool" << endreq;
            return sc;
        }

        if (getCalInfo().isFailure()) {
            log << MSG::ERROR << "Couldn't initialize the CAL constants" << endreq;
            return StatusCode::FAILURE;
        }
    } else {
        return StatusCode::FAILURE;
    }

    // load up the map

    addItem("CalEnergyRaw",  &CAL_EnergyRaw);
    addItem("CalEnergyCorr", &CAL_EnergyCorr);

    addItem("CalLeakCorr",   &CAL_Leak_Corr); 
    addItem("CalEdgeCorr",   &CAL_Edge_Corr);
    addItem("CalTotalCorr",  &CAL_Total_Corr);

    addItem("CalCsIRLn",     &CAL_CsI_RLn);
    addItem("CalTotRLn",     &CAL_Tot_RLn);
    addItem("CalCntRLn",     &CAL_Cntr_RLn);
    addItem("CalLATRLn",     &CAL_LAT_RLn);
    addItem("CalDeadTotRat", &CAL_DeadTot_Rat);
    addItem("CalDeadCntRat", &CAL_DeadCnt_Rat);

    addItem("CalTPred",      &CAL_t_Pred);
    addItem("CalDeltaT",     &CAL_deltaT);

    addItem("CalGapFraction",&CAL_Gap_Fraction);
    addItem("CalTwrEdgeCntr",&CAL_TwrEdgeCntr);
    addItem("CalTwrEdge",    &CAL_TwrEdgeTop);
    addItem("CalLATEdge",    &CAL_LATEdge);
    addItem("CalEdgeEnergy", &CAL_EdgeEnergy);
    addItem("CalTrackDoca",  &CAL_Track_DOCA);
    addItem("CalTrackAngle", &CAL_Track_Angle);
    addItem("CalTrackSep",   &CAL_Track_Sep);

    addItem("CalELayer0",    &CAL_eLayer[0]);
    addItem("CalELayer1",    &CAL_eLayer[1]);
    addItem("CalELayer2",    &CAL_eLayer[2]);
    addItem("CalELayer3",    &CAL_eLayer[3]);
    addItem("CalELayer4",    &CAL_eLayer[4]);
    addItem("CalELayer5",    &CAL_eLayer[5]);
    addItem("CalELayer6",    &CAL_eLayer[6]);
    addItem("CalELayer7",    &CAL_eLayer[7]);
	addItem("CalLghtAsym0",  &CAL_LightAsym[0]);
	addItem("CalLghtAsym1",  &CAL_LightAsym[1]);
	addItem("CalLghtAsym2",  &CAL_LightAsym[2]);
	addItem("CalLghtAsym3",  &CAL_LightAsym[3]);
	addItem("CalLghtAsym4",  &CAL_LightAsym[4]);
	addItem("CalLghtAsym5",  &CAL_LightAsym[5]);
	addItem("CalLghtAsym6",  &CAL_LightAsym[6]);
	addItem("CalLghtAsym7",  &CAL_LightAsym[7]);
    addItem("CalLyr0Ratio",  &CAL_Lyr0_Ratio);
    addItem("CalLyr7Ratio",  &CAL_Lyr7_Ratio);
    addItem("CalBkHalfRatio",&CAL_BkHalf_Ratio);

    addItem("CalXtalsTrunc", &CAL_Num_Xtals_Trunc);
    addItem("CalNumXtals",   &CAL_Num_Xtals);
    addItem("CalXtalRatio",  &CAL_Xtal_Ratio);
    addItem("CalXtalMaxEne", &CAL_Xtal_maxEne);
    addItem("CalMaxNumXtalsInLayer", 
                             &CAL_Max_Num_Xtals_In_Layer);

    addItem("CalTransRms",    &CAL_Trans_Rms);
    addItem("CalLongRms",     &CAL_Long_Rms);
    addItem("CalLRmsAsym",    &CAL_LRms_Asym);

    addItem("CalLongSkew",    &CAL_Long_Skew);
    addItem("CalLongSkewNorm",&CAL_Long_Skew_Norm);

    addItem("CalMIPDiff",   &CAL_MIP_Diff);
    addItem("CalMIPRatio",  &CAL_MIP_Ratio);

    addItem("CalXEcntr",     &CAL_xEcntr);
    addItem("CalYEcntr",     &CAL_yEcntr);
    addItem("CalZEcntr",     &CAL_zEcntr);
    addItem("CalXDir",       &CAL_xdir);
    addItem("CalYDir",       &CAL_ydir);
    addItem("CalZDir",       &CAL_zdir);
    addItem("CalX0",         &CAL_x0);
    addItem("CalY0",         &CAL_y0);
	addItem("CalLyr0X0",      &CAL_X0_Lyr[0]);
	addItem("CalLyr0Y0",      &CAL_Y0_Lyr[0]);
    addItem("CalTopGapDist", &CAL_Top_Gap_Dist);

    addItem("CalXEcntr2",     &CAL_xEcntr2);
    addItem("CalYEcntr2",     &CAL_yEcntr2);
    addItem("CalZEcntr2",     &CAL_zEcntr2);
    addItem("CalXDir2",       &CAL_xdir2);
    addItem("CalYDir2",       &CAL_ydir2);
    addItem("CalZDir2",       &CAL_zdir2);
    addItem("CalPosDirChisq",       &CAL_posdir_chisq);
    addItem("CalPosDirNLayers",       &CAL_posdir_nlayers);
    addItem("CalNSaturated",       &CAL_nsaturated);

    addItem("CalUberEnergy",     &CAL_energy_uber);
    addItem("CalUberXEcntr",     &CAL_xEcntr_uber);
    addItem("CalUberYEcntr",     &CAL_yEcntr_uber);
    addItem("CalUberZEcntr",     &CAL_zEcntr_uber);
    addItem("CalUberXDir",       &CAL_xdir_uber);
    addItem("CalUberYDir",       &CAL_ydir_uber);
    addItem("CalUberZDir",       &CAL_zdir_uber);

    addItem("CalUberXEcntr2",     &CAL_xEcntr2_uber);
    addItem("CalUberYEcntr2",     &CAL_yEcntr2_uber);
    addItem("CalUberZEcntr2",     &CAL_zEcntr2_uber);
    addItem("CalUberXDir2",       &CAL_xdir2_uber);
    addItem("CalUberYDir2",       &CAL_ydir2_uber);
    addItem("CalUberZDir2",       &CAL_zdir2_uber);

    addItem("CalNumClusters",     &CAL_num_clusters);
    addItem("CalRestEnergy",      &CAL_rest_energy);
    addItem("CalRestNumXtals",    &CAL_rest_numXtals);

    addItem("CalTrkXtalRms",       &CAL_track_rms);
    addItem("CalTrkXtalRmsE",      &CAL_track_E_rms);
    addItem("CalTrkXtalRmsTrunc",  &CAL_track_rms_trunc);
    addItem("CalTrkXtalRmsETrunc", &CAL_track_E_rms_trunc);

    addItem("CalCfpEnergy",  &CAL_cfp_energy);
    addItem("CalCfpChiSq",   &CAL_cfp_totChiSq);
    addItem("CalCfpEnergyUB",&CAL_cfp_energyUB);
    addItem("CalCfpEffRLn",  &CAL_cfp_calEffRLn);
    addItem("CalCfpTkrRLn",  &CAL_cfp_tkrRLn);
    addItem("CalCfpAlpha",  &CAL_cfp_alpha);
    addItem("CalCfpTmax",  &CAL_cfp_tmax);
    addItem("CalCfpFitErrFlg",  &CAL_cfp_fiterrflg);
    addItem("CalCfpCalEnergy",  &CAL_cfp_calfit_energy);
    addItem("CalCfpCalChiSq",   &CAL_cfp_calfit_totChiSq);
    addItem("CalCfpCalEffRLn",  &CAL_cfp_calfit_calEffRLn);
    addItem("CalCfpCalAlpha",  &CAL_cfp_calfit_alpha);
    addItem("CalCfpCalTmax",  &CAL_cfp_calfit_tmax);
    addItem("CalCfpCalFitErrFlg",  &CAL_cfp_calfit_fiterrflg);
    //addItem("CalLllEnergy",  &CAL_lll_energy);
    //addItem("CalLllEneErr",  &CAL_lll_energyErr);
    //addItem("CalTklEnergy",  &CAL_tkl_energy);
    //addItem("CalTklEneErr",  &CAL_tkl_energyErr);
    addItem("CalLkHdEnergy", &CAL_LkHd_energy);
    addItem("CalLkHdEneErr", &CAL_LkHd_energyErr);
    addItem("CalLkHdEnergyUB", &CAL_LkHd_energyUB);

    addItem("CalRmsLayerE",  &CAL_rmsLayerE);
    addItem("CalRmsLayerEBack",  &CAL_rmsLayerEBack);
    addItem("CalNLayersRmsBack", &CAL_nLayersRmsBack);
    addItem("CalEAveBack", &CAL_eAveBack);
    addItem("CalLayer0Ratio", &CAL_layer0Ratio);
    addItem("CalXPosRmsLL", &CAL_xPosRmsLastLayer);
    addItem("CalYPosRmsLL", &CAL_yPosRmsLastLayer);

    addItem("CalTrSizeCalT68",&CAL_TS_CAL_T_68);
    addItem("CalTrSizeCalT90",&CAL_TS_CAL_T_90);
    addItem("CalTrSizeCalT95",&CAL_TS_CAL_T_95);
    addItem("CalTrSizeCalT99",&CAL_TS_CAL_T_99);
    addItem("CalTrSizeCalT100",&CAL_TS_CAL_T_100);
    addItem("CalTrSizeCalTL68",&CAL_TS_CAL_TL_68);
    addItem("CalTrSizeCalTL90",&CAL_TS_CAL_TL_90);
    addItem("CalTrSizeCalTL95",&CAL_TS_CAL_TL_95);
    addItem("CalTrSizeCalTL99",&CAL_TS_CAL_TL_99);
    addItem("CalTrSizeCalTL100",&CAL_TS_CAL_TL_100);

    addItem("CalTrSizeCal2T68",&CAL_TS_CAL2_T_68);
    addItem("CalTrSizeCal2T90",&CAL_TS_CAL2_T_90);
    addItem("CalTrSizeCal2T95",&CAL_TS_CAL2_T_95);
    addItem("CalTrSizeCal2T99",&CAL_TS_CAL2_T_99);
    addItem("CalTrSizeCal2T100",&CAL_TS_CAL2_T_100);
//     addItem("CalTrSizeCal2TL68",&CAL_TS_CAL2_TL_68);
//     addItem("CalTrSizeCal2TL90",&CAL_TS_CAL2_TL_90);
//     addItem("CalTrSizeCal2TL95",&CAL_TS_CAL2_TL_95);
//     addItem("CalTrSizeCal2TL99",&CAL_TS_CAL2_TL_99);
//     addItem("CalTrSizeCal2TL100",&CAL_TS_CAL2_TL_100);

    addItem("CalTrSizeTkrT68",&CAL_TS_TKR_T_68);
    addItem("CalTrSizeTkrT90",&CAL_TS_TKR_T_90);
    addItem("CalTrSizeTkrT95",&CAL_TS_TKR_T_95);
    addItem("CalTrSizeTkrT99",&CAL_TS_TKR_T_99);
    addItem("CalTrSizeTkrT100",&CAL_TS_TKR_T_100);
    addItem("CalTrSizeTkrTL68",&CAL_TS_TKR_TL_68);
    addItem("CalTrSizeTkrTL90",&CAL_TS_TKR_TL_90);
    addItem("CalTrSizeTkrTL95",&CAL_TS_TKR_TL_95);
    addItem("CalTrSizeTkrTL99",&CAL_TS_TKR_TL_99);
    addItem("CalTrSizeTkrTL100",&CAL_TS_TKR_TL_100);

    zeroVals();

    //m_ubInterpolate = new UBinterpolate("$(ANALYSISNTUPLEROOT)/calib/BiasMapCalCfpEnergy.txt");
    m_ubInterpolate = new UBinterpolate("$(ANALYSISNTUPLEDATAPATH)/BiasMapCalCfpEnergy.txt");

    return sc;
}

StatusCode CalValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrTrackCol>  
        pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol);
    SmartDataPtr<Event::TkrVertexCol>     
        pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);

    // Recover pointers to CalClusters and Xtals
    SmartDataPtr<Event::CalClusterCol>     
        pCals(m_pEventSvc,EventModel::CalRecon::CalClusterCol);
    SmartDataPtr<Event::CalXtalRecCol> 
        pxtalrecs(m_pEventSvc,EventModel::CalRecon::CalXtalRecCol);

    // Recover pointer to CalEventEnergy info  
#ifdef PRE_CALMOD
    Event::CalEventEnergy* calEventEnergy = 
        SmartDataPtr<Event::CalEventEnergy>(m_pEventSvc, EventModel::CalRecon::CalEventEnergy);
#else
    Event::CalEventEnergyCol * calEventEnergyCol = 
        SmartDataPtr<Event::CalEventEnergyCol>(m_pEventSvc,EventModel::CalRecon::CalEventEnergyCol) ;
    Event::CalEventEnergy * calEventEnergy = 0 ;
    if ((calEventEnergyCol!=0)&&(!calEventEnergyCol->empty()))
        calEventEnergy = calEventEnergyCol->front() ;
#endif


    // If calEventEnergy then fill TkrEventParams
    // Note: TkrEventParams initializes to zero in the event of no CalEventEnergy
    double m_radLen_Stuff     = 0.;
    double m_radLen_CntrStuff = 0.; 
    if (calEventEnergy != 0)
    {
        // Extraction of results from CalValCorrTool in CalRecon... 
        for(Event::CalCorToolResultCol::iterator corIter = calEventEnergy->begin(); corIter != calEventEnergy->end(); corIter++)
        {
            Event::CalCorToolResult corResult = **corIter;

            if (corResult.getCorrectionName() == "CalValsCorrTool")
            {
                CAL_EnergyCorr   = corResult["CorrectedEnergy"];
                CAL_Leak_Corr    = corResult["LeakCorrection"];
                CAL_Edge_Corr    = corResult["EdgeCorrection"]; 
                CAL_Gap_Fraction = corResult["GapFraction"];     
                CAL_Total_Corr  = corResult["TotalCorrection"]; 

                CAL_CsI_RLn        = corResult["CsIRLn"];
                CAL_LAT_RLn        = corResult["LATRLn"]; 
                m_radLen_Stuff     = corResult["StuffRLn"];
                m_radLen_CntrStuff = corResult["CntrRLnStuff"];
                CAL_Cntr_RLn       = corResult["CntrRLn"];
                CAL_t_Pred  = corResult["PredCntr"]; 
                CAL_deltaT  = corResult["DeltaCntr"];
                CAL_Tot_RLn = corResult["CALRLn"];

               // CAL_x0 = corResult["CalTopX0"]; See below under Imaging Calorimeter Top
               // CAL_y0 = corResult["CalTopY0"];

                //CAL_RmsE = corResult["RmsE"];
                //CAL_numLayersRms = corResult["NumLayersRms"];
            }
            else if (corResult.getCorrectionName() == "CalFullProfileTool")
            {
                CAL_cfp_energy    = corResult.getParams().getEnergy();
                CAL_cfp_totChiSq  = corResult["totchisq"];
                CAL_cfp_calEffRLn = corResult["cal_eff_RLn"];
                CAL_cfp_tkrRLn = corResult["tkr_RLn"];
                CAL_cfp_alpha = corResult["alpha"];
                CAL_cfp_tmax = corResult["tmax"];
                CAL_cfp_fiterrflg = corResult["fitflag"];
                CAL_cfp_calfit_energy    = corResult["calfit_fit_energy"];
                CAL_cfp_calfit_totChiSq  = corResult["calfit_totchisq"];
                CAL_cfp_calfit_calEffRLn = corResult["calfit_cal_eff_RLn"];
                CAL_cfp_calfit_alpha = corResult["calfit_alpha"];
                CAL_cfp_calfit_tmax = corResult["calfit_tmax"];
                CAL_cfp_calfit_fiterrflg = corResult["calfit_fitflag"];
                if ( CAL_cfp_energy == 0 ) {
                    CAL_cfp_energyUB = 0;
                    //                    std::cout << "CalValsTool CAL_cfp_energy " << CAL_cfp_energy << ' ' << CAL_cfp_energyUB << std::endl;
                }
                else {
                    float tkr1ZDir = -1;
                    if ( m_pTkrTool->getVal("Tkr1ZDir", tkr1ZDir).isSuccess() ) {
                        if ( tkr1ZDir == 0 )
                            tkr1ZDir = -1;
                    }
                    const float bias = m_ubInterpolate->interpolate(log10(CAL_cfp_energy), tkr1ZDir);
                    CAL_cfp_energyUB = bias == 0 ? -1 : CAL_cfp_energy / bias;
                    //                    std::cout << "EvtValsTool CAL_cfp_energy " << CAL_cfp_energy << " ( " << log10(CAL_cfp_energy) << " ) " << CAL_cfp_energyUB << " ( " << log10(CAL_cfp_energyUB) << " ) " << tkr1ZDir << std::endl;
                }
            }
            // Removed 5/5/09 LSR
            //else if (corResult.getCorrectionName() == "CalLastLayerLikelihoodTool")
            //{
            //    CAL_lll_energy    = corResult.getParams().getEnergy();
            //    CAL_lll_energyErr = corResult.getParams().getEnergyErr();
            //}
            //else if (corResult.getCorrectionName() == "CalTkrLikelihoodTool")
            //{
            //    CAL_tkl_energy    = corResult.getParams().getEnergy();
            //    CAL_tkl_energyErr = corResult.getParams().getEnergyErr();
            //}
            else if (corResult.getCorrectionName() == "CalLikelihoodManagerTool")
            {
                CAL_LkHd_energy    = corResult.getParams().getEnergy();
                CAL_LkHd_energyErr = corResult.getParams().getEnergyErr();
                CAL_LkHd_energyUB = CAL_LkHd_energy / ( 1.003 - 0.005345 * log10(CAL_LkHd_energy) );
                //                std::cout << "CalValsTool CAL_LkHd_energy " << CAL_LkHd_energy << ' ' << CAL_LkHd_energyUB << std::endl;
            }
        }
    }

    //Make sure we have valid cluster data and some energy
    if (!pCals) return sc;
    if (pCals->empty()) return sc;

    // Extract the uber information which is located at the end of the list
    Event::CalCluster* calCluster = pCals->back();

    CAL_energy_uber = calCluster->getCalParams().getEnergy();

    CAL_xEcntr_uber   = calCluster->getPosition().x();
    CAL_yEcntr_uber   = calCluster->getPosition().y();
    CAL_zEcntr_uber   = calCluster->getPosition().z();
    CAL_xdir_uber     = calCluster->getDirection().x();
    CAL_ydir_uber     = calCluster->getDirection().y();
    CAL_zdir_uber     = calCluster->getDirection().z();

    // Get pos and dir determined using only the transverse position information
    CAL_xEcntr2_uber  = calCluster->getFitParams().getCentroid().x();
    CAL_yEcntr2_uber  = calCluster->getFitParams().getCentroid().y();
    CAL_zEcntr2_uber  = calCluster->getFitParams().getCentroid().z();
    CAL_nsaturated    = calCluster->getNumSaturatedXtals();
    CAL_xdir2_uber    = calCluster->getFitParams().getAxis().x();
    CAL_ydir2_uber    = calCluster->getFitParams().getAxis().y();
    CAL_zdir2_uber    = calCluster->getFitParams().getAxis().z();

    CAL_num_clusters  = pCals->size();
    CAL_rest_energy   = 0.;
    CAL_rest_numXtals = 0.;

    // If the collection contains more than one entry then we have multiple
    // clusters. If we don't then our current pointer points to the lone cluster
    if (pCals->size() > 1)
    {
        Event::CalClusterCol::iterator calClusIter = pCals->begin();
        calCluster = *calClusIter++;

        // Loop through remaining crystals and add up the energy
        while(calClusIter != pCals->end())
        {
            Event::CalCluster* cluster = *calClusIter++;

            CAL_rest_energy += cluster->getCalParams().getEnergy();
        }

        // Remove the uber cluster energy
        CAL_rest_energy -= CAL_energy_uber;
    }

    CAL_EnergyRaw  = calCluster->getCalParams().getEnergy();
    if(CAL_EnergyRaw<1.0) return sc;

    for(int i = 0; i<m_nLayers; i++) CAL_eLayer[i] = (*calCluster)[i].getEnergy();

    CAL_Trans_Rms = sqrt(calCluster->getRmsTrans()/CAL_EnergyRaw);

    float logRLn = 0; 
    if ((CAL_LAT_RLn - CAL_Cntr_RLn) > 0.0) {
        logRLn = log(CAL_LAT_RLn - CAL_Cntr_RLn);
    }
    if (logRLn > 0.0) {
        CAL_Long_Rms  = sqrt(calCluster->getRmsLong()/CAL_EnergyRaw) / logRLn;
    }
    CAL_LRms_Asym = calCluster->getRmsLongAsym();

    // Start with the skewness stuff.
    // At this poinu we need to normalize to the sum of weights and divide by sigma**3.
    CAL_Long_Skew = calCluster->getSkewnessLong();
    CAL_Long_Skew_Norm = -9999.;
    if (CAL_EnergyRaw > 0.0) {
      double sigma = sqrt(calCluster->getRmsLong()/CAL_EnergyRaw);
      if (sigma > 0.0){
	// First normalize for the sum of weights and to the third power of the longitudinal RMS.
    	CAL_Long_Skew = CAL_Long_Skew/(CAL_EnergyRaw*sigma*sigma*sigma);
	// Then calculate what we would expect for an EM shower (downward going).
	// This involves a numerical integral, with the basic formulas from the PDG.
	int numSteps = 6;                             // Number of steps for the trapezium integral.
	double Ec = 11.17;                            // Critical energy for CsI (in MeV).
	double tmin = CAL_LAT_RLn - CAL_CsI_RLn; 
	double tmax = CAL_LAT_RLn;
	double b  = 0.5;                              // b = 0.5 seems appropriate, here.
	double c  = -0.5;                             // c = -0.5 for electrons, c = 0.5 for gammas.
	double a  = b*(log(CAL_EnergyCorr/Ec) - 0.5); // Use CalEnergyCorr. Can we improve?
	double t1 = tmin;
	double t2 = 0.0;
	double p1 = 0.0;
	double p2 = 0.0;
	double norm = 0.0;
	double mom1 = 0.0;
	double mom2 = 0.0;
	double mom3 = 0.0;
	double stepSize = (tmax - tmin)/numSteps;
	for (int step = 0; step < numSteps; step++) {
	  t2 = tmin + (step + 1)*stepSize;
	  p1 = pow(t1, a) * exp(-b*t1);
	  p2 = pow(t2, a) * exp(-b*t2);
	  norm += 0.5 * (p1 + p2);
	  mom1 += 0.5 * (p1*t1 + p2*t2);
	  mom2 += 0.5 * (p1*t1*t1 + p2*t2*t2);
	  mom3 += 0.5 * (p1*t1*t1*t1 + p2*t2*t2*t2);
	  t1 = t2;
	}
	norm *= stepSize;
	mom1 *= stepSize;
	mom2 *= stepSize;
	mom3 *= stepSize;
	mom1 /= norm;
	mom2 /= norm;
	mom3 /= norm;
	mom2 = mom2 - mom1*mom1;
	double gamma = (mom3 - 3*mom1*mom2 - mom1*mom1*mom1)/pow(mom2, 1.5);
	if (gamma != 0.0) CAL_Long_Skew_Norm = -CAL_Long_Skew/gamma;
      }
    }

    if(CAL_EnergyRaw>0.0) {
        CAL_Lyr0_Ratio  = CAL_eLayer[0]/CAL_EnergyRaw;
        if(m_nLayers==8) {
            CAL_Lyr7_Ratio  = CAL_eLayer[7]/CAL_EnergyRaw;
            CAL_BkHalf_Ratio = (CAL_eLayer[4]+CAL_eLayer[5]+
                CAL_eLayer[6]+CAL_eLayer[7])/CAL_EnergyRaw;
        }
    }

    int no_xtals=0;
    CAL_Xtal_maxEne = 0.;

    // Local array for #xtals in layer with the most xtals (cut on e>=5MeV);
    std::vector<int> xtalCount(m_nLayers,0);
	std::vector<double> lightAsym(m_nLayers,0.);
	std::vector<double> eneLogLayer(m_nLayers,0.);
	std::vector<double> logX0(m_nLayers,0.);
	std::vector<double> logY0(m_nLayers, 0.);

    Event::CalXtalRecCol::const_iterator jlog;
    if (pxtalrecs) {
        // Find Xtal with max. energy
        CAL_Num_Xtals = (float)pxtalrecs->size();
        for( jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog) {
            const Event::CalXtalRecData& recLog = **jlog;    
            double eneLog = recLog.getEnergy();
			double enePos = recLog.getEnergy(0, idents::CalXtalId::POS);
			double eneNeg = recLog.getEnergy(0, idents::CalXtalId::NEG);
			Point pos = recLog.getPosition();
            if(enePos > CAL_Xtal_maxEne) CAL_Xtal_maxEne = enePos;
			if(eneNeg > CAL_Xtal_maxEne) CAL_Xtal_maxEne = eneNeg;
            idents::CalXtalId xtalId = recLog.getPackedId();
            int layer = xtalId.getLayer();
            if(eneLog>5.0) xtalCount[layer]++;

			// Light asymetry section for testing

			if(eneLog > 5 && eneLog > eneLogLayer[layer]) {
				eneLogLayer[layer] = eneLog;
				lightAsym[layer] = (enePos - eneNeg)/(enePos + eneNeg);	
				logX0[layer] = pos.x();
				logY0[layer] = pos.y();
			}
        }
		for(int il=0; il<8; il++) {
			CAL_LightAsym[il]=lightAsym[il];
			CAL_X0_Lyr[il] = logX0[il];
			CAL_Y0_Lyr[il] = logY0[il];
		}

        std::vector<int>::const_iterator itC = 
            std::max_element(xtalCount.begin(),xtalCount.end());
        CAL_Max_Num_Xtals_In_Layer = *itC;

        // Number of Xtals
        no_xtals=pxtalrecs->size();
    }
    int no_xtals_trunc=calCluster->getNumTruncXtals();
    CAL_Xtal_Ratio= (no_xtals>0) ? float(no_xtals_trunc)/no_xtals : 0;
    CAL_Num_Xtals_Trunc = float(no_xtals_trunc); 

    // No use in continuing if too little energy in CAL
    if(CAL_EnergyRaw < 5.) return sc;  

    Point  cal_pos  = calCluster->getPosition();
    Vector cal_dir  = calCluster->getDirection();

    CAL_xEcntr      = cal_pos.x();
    CAL_yEcntr      = cal_pos.y();
    CAL_zEcntr      = cal_pos.z();
    CAL_xdir        = cal_dir.x();
    CAL_ydir        = cal_dir.y();
    CAL_zdir        = cal_dir.z();

    // Get pos and dir determined using only the transverse position information
    CAL_xEcntr2        = calCluster->getFitParams().getCentroid().x();
    CAL_yEcntr2        = calCluster->getFitParams().getCentroid().y();
    CAL_zEcntr2        = calCluster->getFitParams().getCentroid().z();
    CAL_posdir_chisq   = calCluster->getFitParams().getChiSquare();
    CAL_posdir_nlayers = calCluster->getFitParams().getFitLayers();
    CAL_nsaturated     = calCluster->getNumSaturatedXtals();
    CAL_xdir2          = calCluster->getFitParams().getAxis().x();
    CAL_ydir2          = calCluster->getFitParams().getAxis().y();
    CAL_zdir2          = calCluster->getFitParams().getAxis().z();

    // Get the lower and upper limits for the CAL in the installed towers
    double deltaX = 0.5*(m_xNum*m_towerPitch - m_calXWidth);
    double deltaY = 0.5*(m_yNum*m_towerPitch - m_calYWidth);
    double calXLo = m_tkrGeom->getLATLimit(0,LOW)  + deltaX;
    double calXHi = m_tkrGeom->getLATLimit(0,HIGH) - deltaX;
    double calYLo = m_tkrGeom->getLATLimit(1,LOW)  + deltaY;
    double calYHi = m_tkrGeom->getLATLimit(1,HIGH) - deltaY;

    
    // collect the CAL edge energy
    // the edge is larger of the width and length of the layer
    // we can do better if this is at all useful.
    if (pxtalrecs) {

        // do the Cal[X/Y]PosRmsLL here
        double eRmsLL = 0;
        double xLL = 0;
        double xSqLL = 0;
        double yLL = 0;
        double ySqLL = 0;
        int    nRmsLL = 0;
        bool   doLL = (m_nLayers>1 && m_nCsI>1); // true for Glast, false for EGRET

        for( jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog) {
            const Event::CalXtalRecData& recLog = **jlog;    
            Point pos = recLog.getPosition();
                double eneLog = recLog.getEnergy();
            double xPos = pos.x(); double yPos = pos.y();
            double minX, minY;
            minX = std::min(fabs(xPos-calXLo),fabs(xPos-calXHi));
            minY = std::min(fabs(yPos-calYLo),fabs(yPos-calYHi));

            // for last-layer rms
            idents::CalXtalId id = recLog.getPackedId();
            int layer = id.getLayer();
            if (eneLog>0 && layer==m_nLayers-1 && doLL) {
                nRmsLL++;
                eRmsLL += eneLog;
                xLL  += eneLog*xPos;
                yLL  += eneLog*yPos;
                xSqLL += eneLog*xPos*xPos;
                ySqLL += eneLog*yPos*yPos;
            }

            // for later... no time to check this out now!
            /*
            idents::CalXtalId id = recLog.getPackedId();
            if (id.isX()) {
            double minX = 
            std::min(fabs(xPos-(calXLo+m_deltaSide)),fabs(xPos-(calXHi-m_deltaSide)));
            double minY = std::min(fabs(yPos-calYLo),fabs(yPos-calYHi));
            } else {
            double minX = std::min(fabs(xPos-calXLo),fabs(xPos-calXHi));
            double minY = 
            std::min(fabs(yPos-(calYLo+m_deltaSide)),fabs(yPos-(calYHi-m_deltaSide)));
            }
            */
            if ((minX<_deltaEdge) ||  (minY<_deltaEdge)) {
                double eneLog = recLog.getEnergy();
                CAL_EdgeEnergy += eneLog;
            }
        }

        // finish last-layer rms 
        if(nRmsLL>1) {
            double xAveLL = xLL/eRmsLL;
            double xVarLL = std::max(0.0, xSqLL/eRmsLL - xAveLL*xAveLL);
            double yAveLL = yLL/eRmsLL;
            double yVarLL = std::max(0.0, ySqLL/eRmsLL - yAveLL*yAveLL);
            CAL_xPosRmsLastLayer = (float) sqrt( xVarLL*nRmsLL/(nRmsLL-1));
            CAL_yPosRmsLastLayer = (float) sqrt( yVarLL*nRmsLL/(nRmsLL-1));
        }
    }

    // Compare Track (or vertex) to Cal soltion
    Point x0  =  cal_pos;
    Vector t0 = -cal_dir;
    if(calCluster->getRmsLong() < .1 ) { // Trap no-calc. condition
        x0 = Point(0., 0., 0.);
        t0 = Vector(0., 0., -1.);
    }

    // If there are tracks, use the best one
    int num_tracks = 0;
    Point x1, x2; 
    Vector t1, t2; 

    if(pTracks) num_tracks = pTracks->size();

    Event::TkrTrackColConPtr pTrack1;
    Event::TkrTrack* track_1 = 0;

    if(num_tracks > 0) { 
        // Get the first track
        pTrack1 = pTracks->begin();
        track_1 = *pTrack1;

        // Get the start and direction 
        x1 = track_1->getInitialPosition();
        t1 = track_1->getInitialDirection();
        x0 = x1; 
        t0 = t1;

        // Find the distance for energy centroid to track axis
        /* old calculation
        Vector x_diff = x0 - cal_pos;
        double x_diff_sq = x_diff*x_diff;
        double x_diff_t0 = x_diff*t0;
        CAL_Track_DOCA = sqrt(std::max(0.0, x_diff_sq - x_diff_t0*x_diff_t0));
        */
        Doca track1(x1, t1);
        CAL_Track_DOCA = (float)track1.docaOfPoint(cal_pos);

		// Image the top of the Calorimeter - use last hit FILTERED
		const Event::TkrTrackParams& tkr1_params = track_1->back()->getTrackParams(Event::TkrTrackHit::FILTERED);
		double xSlp = tkr1_params.getxSlope();
		double ySlp = tkr1_params.getySlope();
		Vector LastDir = Vector(-xSlp, -ySlp, -1).unit();
		Point LastHit = track_1->back()->getPoint(Event::TkrTrackHit::FILTERED);
		
		//double calZTop = -46.;
		double deltaZ  = m_calZTop - LastHit.z();
		double topArcLength = deltaZ/LastDir.z();
		Point calTopLoc = LastHit + topArcLength*LastDir; 
		Point calMidLayer0Loc = LastHit + (topArcLength - 9.5/LastDir.z())*LastDir; 
		CAL_x0 = calMidLayer0Loc.x();
		CAL_y0 = calMidLayer0Loc.y(); 

		// Event axis locations relative to towers and LAT
    Point pos0(CAL_x0, CAL_y0, m_calZTop);
    int view;
    CAL_TwrEdgeTop  = activeDist(pos0, view);
    CAL_TwrEdgeCntr = activeDist(cal_pos, view);
    // Find the distance closest to an edge
    double dX = std::max(calXLo-CAL_x0, CAL_x0-calXHi);
    double dY = std::max(calYLo-CAL_y0, CAL_y0-calYHi);
    CAL_LATEdge = -std::max(dX, dY);


		// Now create an active distance type variable using the tower pitch
		double integer_part;
		double deltaX_crack = modf(fabs(CAL_x0)/m_towerPitch, &integer_part); 
		if(deltaX_crack > .5) deltaX_crack = 1. - deltaX_crack;
		double deltaY_crack = modf(fabs(CAL_y0)/m_towerPitch, &integer_part);
		if(deltaY_crack > .5) deltaY_crack = 1. - deltaY_crack;
		CAL_Top_Gap_Dist = m_towerPitch*std::min(deltaX_crack, deltaY_crack);

        if(num_tracks > 1) { 
            // Get the second track
            pTrack1++;
            Event::TkrTrack* track_1 = *pTrack1;

            // Get the start and direction 
            x2 = track_1->getInitialPosition();
            t2 = track_1->getInitialDirection();

            Point x2Top = x2 + (pos0.z()-x2.z())/t2.z() * t2; 
            CAL_Track_Sep = (x2Top - pos0).mag();
        }


        // If vertexed - use first vertex
        if(pVerts && pVerts->size()>0) {
            Event::TkrVertex* vertex = *(pVerts->begin()); 
            x0 = vertex->getPosition();
            t0 = vertex->getDirection();
        }

    // try Bill's new vars... 
        if (pxtalrecs) {
            //make a map of xtal energy by doca
            // we need this so that we can drop the largest entries for the truncated calc.
            std::multimap<double, double> docaMap;
            std::multimap<double, double>::iterator mapIter;

            for( jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog) {
                const Event::CalXtalRecData& recLog = **jlog;    
                Point pos = recLog.getPosition();
                double doca = track1.docaOfPoint(pos);
                double energy = recLog.getEnergy();

                std::pair<double, double> docaPair(doca, energy);
                docaMap.insert(docaPair); 
            }

            unsigned int mapSize = docaMap.size();
            int truncSize = std::min((int)((mapSize*truncFraction)+0.5),(int) mapSize);

            // make up the vars by iterating over the map
            // for the truncated var, stop after truncFraction of the xtals
            int mapCount = 0;
            // 4 measures of the Track-Cal RMS
            float totalE1 = 0.0;
            float totalE2 = 0.0;
            for(mapIter=docaMap.begin(); mapIter!=docaMap.end();++mapIter) {
                ++mapCount;
                double doca = mapIter->first;
                double ene  = mapIter->second;
                double docaSq = doca*doca;
                CAL_track_rms += (float) docaSq;
                CAL_track_E_rms += (float) docaSq*ene;
                totalE1 += (float) ene;
                if (mapCount<=truncSize) {
                    CAL_track_rms_trunc += (float) docaSq;
                    CAL_track_E_rms_trunc += (float) docaSq*ene;
                    totalE2 += (float) ene;
                }     
            }
            CAL_track_rms = sqrt(CAL_track_rms/mapSize);
            if (totalE1>0.0) {
                CAL_track_E_rms = 
                    sqrt(std::max(0.0f, CAL_track_E_rms)/totalE1);
            }
            CAL_track_rms_trunc = sqrt(CAL_track_rms_trunc/truncSize);
            if (totalE2>0.0) {
                CAL_track_E_rms_trunc = 
                    sqrt(std::max(0.0f, CAL_track_E_rms_trunc)/totalE2);
            }
        }
    }

    // Find angle between Track and Cal. Moment axis
    // Note: the direction in Cal is opposite to tracking!
    if(num_tracks>0 && fabs(cal_dir.x()) < 1.) {
        double cosCalt0 = std::min(1., -t0*cal_dir); 
        cosCalt0 = std::max(cosCalt0, -1.);  // just in case...
        CAL_Track_Angle = acos(cosCalt0);
    } else {
        CAL_Track_Angle = -.1f; 
    }

    // Energy compared to muon 
    CAL_MIP_Diff    = CAL_EnergyRaw- 12.07*CAL_CsI_RLn;
    const double minRadLen = 0.1;
    CAL_MIP_Ratio   = CAL_EnergyRaw/(12.07*std::max(CAL_CsI_RLn*1., minRadLen));

    // Ratios of rad. len. in material other then CsI
    CAL_DeadTot_Rat = m_radLen_Stuff/std::max(minRadLen, CAL_Tot_RLn*1.);
    CAL_DeadCnt_Rat = m_radLen_CntrStuff/std::max(minRadLen, CAL_Cntr_RLn*1.); 

    // Here we do the CAL_RmsE calculation

    Point xEnd;
    Vector tEnd;
    double arclen;

    try {
    // get the last point on the best track
    if(num_tracks>0) {
        Event::TkrTrackHitVecConItr hitIter = (*track_1).end();
        hitIter--;
        const Event::TkrTrackHit* hit = *hitIter;
        xEnd = hit->getPoint(Event::TkrTrackHit::SMOOTHED);
        tEnd = hit->getDirection(Event::TkrTrackHit::SMOOTHED);
        double eCosTheta = fabs(tEnd.z());

        // find the parameters to start the swim
        double deltaZ = xEnd.z() - m_calZBot;
        arclen   = fabs(deltaZ/eCosTheta);

        // do the swim
        m_G4PropTool->setStepStart(xEnd, tEnd);
        m_G4PropTool->step(arclen);
    
        // collect the radlens by layer
        int numSteps = m_G4PropTool->getNumberSteps();
        std::vector<double> rlCsI(m_nLayers, 0.0);
        std::vector<bool>   useLayer(m_nLayers, true);
        idents::VolumeIdentifier volId;
        idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();
        int istep  = 0;
        for(; istep < numSteps; ++istep) {
            volId = m_G4PropTool->getStepVolumeId(istep);
            volId.prepend(prefix);
            bool inXtal = ( volId.size()>7 && volId[0]==0 
                && volId[3]==0 && volId[7]==0 ? true : false );
            if(inXtal) {
                int layer = volId[4];
                double radLen_step = m_G4PropTool->getStepRadLength(istep);
                rlCsI[layer] += radLen_step;
            }
        }

        double eTot = 0;
        double eNorm0 = 0;
        double e2   = 0;
        int CAL_nLayersRms = 0;
        int layer;

        for (layer=0; layer<m_nLayers; ++layer) {
            if (rlCsI[layer]<0.5) useLayer[layer] = false;
            if (CAL_eLayer[layer]<0.05*CAL_EnergyRaw || CAL_EnergyRaw<=0) {
                useLayer[layer] = false;
            }
        }

        for (layer=0; layer<m_nLayers; ++layer) {
            if(!useLayer[layer]) continue;
            CAL_nLayersRms++;
            double eNorm = CAL_eLayer[layer]/rlCsI[layer];
            if(layer==0) { 
                eNorm0 = eNorm;
            } else {
                CAL_nLayersRmsBack++;
            }
            eTot += eNorm;
            e2 +=   eNorm*eNorm;
            //std::cout << "layer " << layer << ", r.l. " << rlCsI[layer] 
            //    << ", eLayer " << CAL_eLayer[layer] 
            //    << ", rlNorm " << rlCsI[layer]*eCosTheta) << ", eL_norm " 
            //    << eNorm << std::endl;
        }

        double eAve = 0;
        double eAveBack = 0;
        if(CAL_nLayersRms>1) {
            eAve = eTot/CAL_nLayersRms;
            CAL_rmsLayerE = 
                (float)sqrt((e2 - CAL_nLayersRms*eAve*eAve)/(CAL_nLayersRms-1))/eAve;
            //std::cout << "eRms " << CAL_rmsLayerE/eAve 
            //    << ", " << CAL_nLayersRms << " layers" << std::endl;
        }
        if(CAL_nLayersRmsBack>1) {
            eAveBack = (eTot-eNorm0)/(CAL_nLayersRmsBack);
            CAL_rmsLayerEBack = 
                (float) sqrt((e2 - eNorm0*eNorm0 - (CAL_nLayersRmsBack)*eAveBack*eAveBack)
                /(CAL_nLayersRmsBack-1))/eAve;
            //std::cout << "eRmsTrunc " << CAL_rmsLayerEBack/eAveBack 
            //    << ", " << nLayersRmsBack << " layers"<<std::endl;
            CAL_eAveBack = eAveBack;
            CAL_layer0Ratio = eNorm0/eAveBack;
        }
    }
    } catch( std::exception& /*e*/ ) {
        MsgStream log(msgSvc(), name());
        printHeader(log);
        setAnaTupBit();
        log << "See previous exception message." << endreq;
        log << " Skipping the Cal_rmsE calculation and resetting" << endreq;

        CAL_layer0Ratio = _badFloat;
        CAL_eAveBack    = _badFloat;
        CAL_nLayersRmsBack = _badInt;
        CAL_rmsLayerEBack = _badFloat;
    } catch(...) {
        MsgStream log(msgSvc(), name());
        printHeader(log);
        setAnaTupBit();
        log << "Unknown exception: see previous exception message, if any." << endreq;
        log << " Skipping the Cal_rmsE calculation" << endreq;
        log << "Initial track parameters: pos: " << xEnd << endreq 
            << "dir: " << tEnd << " arclen: " << arclen << endreq;

        CAL_layer0Ratio = _badFloat;
        CAL_eAveBack    = _badFloat;
        CAL_nLayersRmsBack = _badInt;
        CAL_rmsLayerEBack = _badFloat;
    }

    //
    // Estimations of the transverse size (Philippe Bruel)
    //
    
    //
    // Perform the estimation with main axis = cal axis
    //
    TSaxisP = calCluster->getPosition();
    Vector TSaxisVin = calCluster->getDirection();
    TSaxisV = TSaxisVin.unit();
    //
    // Filling TSdist...
    //
    TSfillTSdist(pxtalrecs);

    //int i;

    if(TSnlog>0)
      {
	//
	// fill CAL_TS_CAL_T_
	//
	TSfillTS(0);
	CAL_TS_CAL_T_68 = (float)TSgetinterpolationTS(0.68);
	CAL_TS_CAL_T_90 = (float)TSgetinterpolationTS(0.90);
	CAL_TS_CAL_T_95 = (float)TSgetinterpolationTS(0.95);
	CAL_TS_CAL_T_99 = (float)TSgetinterpolationTS(0.99);
	CAL_TS_CAL_T_100 = (float)TSTS[TSnlog-1];
	//
	// fill CAL_TS_CAL_TL_
	//
	TSfillTS(1);
	CAL_TS_CAL_TL_68 = (float)TSgetinterpolationTS(0.68);
	CAL_TS_CAL_TL_90 = (float)TSgetinterpolationTS(0.90);
	CAL_TS_CAL_TL_95 = (float)TSgetinterpolationTS(0.95);
	CAL_TS_CAL_TL_99 = (float)TSgetinterpolationTS(0.99);
	CAL_TS_CAL_TL_100 = (float)TSTS[TSnlog-1];
      }

    //
    // Perform the estimation with centroid and axis determined only with transverse information
    //
    if(CAL_posdir_nlayers>=4)
      {
	TSaxisP = Point(CAL_xEcntr2,CAL_yEcntr2,CAL_zEcntr2);
	TSaxisVin = Vector(CAL_xdir2,CAL_ydir2,CAL_zdir2);
	TSaxisV = TSaxisVin.unit();
	//
	// Filling TSdist...
	//
	TSfillTSdist(pxtalrecs);
	
	if(TSnlog>0)
	  {
	    //
	    // fill CAL_TS_CAL2_T_
	    //
	    TSfillTS(0);
	    CAL_TS_CAL2_T_68 = (float)TSgetinterpolationTS(0.68);
	    CAL_TS_CAL2_T_90 = (float)TSgetinterpolationTS(0.90);
	    CAL_TS_CAL2_T_95 = (float)TSgetinterpolationTS(0.95);
	    CAL_TS_CAL2_T_99 = (float)TSgetinterpolationTS(0.99);
	    CAL_TS_CAL2_T_100 = (float)TSTS[TSnlog-1];
	    //
	    // fill CAL_TS_CAL2_TL_
	    //
	    TSfillTS(1);
	    CAL_TS_CAL2_TL_68 = (float)TSgetinterpolationTS(0.68);
	    CAL_TS_CAL2_TL_90 = (float)TSgetinterpolationTS(0.90);
	    CAL_TS_CAL2_TL_95 = (float)TSgetinterpolationTS(0.95);
	    CAL_TS_CAL2_TL_99 = (float)TSgetinterpolationTS(0.99);
	    CAL_TS_CAL2_TL_100 = (float)TSTS[TSnlog-1];
	  }
      }

    //
    // Perform the estimation with main axis = best track
    //
    int mynumtracks = 0;
    if(pTracks) mynumtracks = pTracks->size();
    if(mynumtracks>0)
      { 
        // Get the first track
        pTrack1 = pTracks->begin();
        track_1 = *pTrack1;
        // Get the start and direction 
        TSaxisP = track_1->getInitialPosition();
        TSaxisVin = track_1->getInitialDirection();
	TSaxisV = TSaxisVin.unit();
	//
	// Filling TSdist...
	//
	TSfillTSdist(pxtalrecs);

	if(TSnlog>0)
	  {
	    //
	    // fill CAL_TS_TKR_T_
	    //
	    TSfillTS(0);
	    CAL_TS_TKR_T_68 = (float)TSgetinterpolationTS(0.68);
	    CAL_TS_TKR_T_90 = (float)TSgetinterpolationTS(0.90);
	    CAL_TS_TKR_T_95 = (float)TSgetinterpolationTS(0.95);
	    CAL_TS_TKR_T_99 = (float)TSgetinterpolationTS(0.99);
	    CAL_TS_TKR_T_100 = (float)TSTS[TSnlog-1];
	    //
	    // fill CAL_TS_TKR_TL_
	    //
	    TSfillTS(1);
	    CAL_TS_TKR_TL_68 = (float)TSgetinterpolationTS(0.68);
	    CAL_TS_TKR_TL_90 = (float)TSgetinterpolationTS(0.90);
	    CAL_TS_TKR_TL_95 = (float)TSgetinterpolationTS(0.95);
	    CAL_TS_TKR_TL_99 = (float)TSgetinterpolationTS(0.99);
	    CAL_TS_TKR_TL_100 = (float)TSTS[TSnlog-1];
	  }
      }

      // throw an exception
      //int ii = 1;
      //int j = 0;
      //int k = ii/j;
      //k++;
    
    return sc;
}

StatusCode CalValsTool::getCalInfo()
{
    m_calZTop = m_tkrGeom->calZTop();
    m_calZBot = m_tkrGeom->calZBot();
    m_calXWidth = m_tkrGeom->calXWidth();
    m_calYWidth = m_tkrGeom->calYWidth();
    /*
    double calModuleWid, calXtalLen;
    m_detSvc->getNumericConstByName("CalModuleWidth", &calModuleWid);
    m_detSvc->getNumericConstByName("CsILength", &calXtalLen);
    double m_deltaSide = calModuleWid - calXtalLen;
    */

    return StatusCode::SUCCESS;
}

double CalValsTool::activeDist(Point pos, int &view) const
{
    double edge = 0.;
    double x = pos.x();
    double y = pos.y();
    double x_twr = globalToLocal(x, m_towerPitch, m_xNum);
    double y_twr = globalToLocal(y, m_towerPitch, m_yNum);

    if( fabs(x_twr) > fabs(y_twr) ) {
        edge = m_towerPitch/2. - fabs(x_twr);
        view = 0; 
    }
    else {
        edge = m_towerPitch/2. - fabs(y_twr);
        view = 1;
    }
    return edge;
}

int CalValsTool::TSfillTSdist(Event::CalXtalRecCol *pxtalrecs)
{
  int i;
  TSnlog = 0;
  for(i=0;i<1536;++i)
    {
      TSdistTL[i] = 0;
      TSdistT[i] = 0;
      TSenergy[i] = 0;
    }

  if(pxtalrecs==NULL) return 1;

  int itow,ilay,icol,itowx,itowy;
  double lambda;
  Point TSxtalP;
  Point TSxtalC;
  Vector TSxtalV;
  Vector TSTC;
  double lambdamax = 326./2;
 
  Event::CalXtalRecCol::const_iterator jlog;
 
  for( jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog)
    {
      const Event::CalXtalRecData& recLog = **jlog;    
      TSxtalP = recLog.getPosition();
      TSenergy[TSnlog] = recLog.getEnergy();
      //
      // computing TSdistTL : using both the transverse and longitudinal position measurements
      //
      TSTC = TSxtalP-TSaxisP;
      TSdistTL[TSnlog] = TSTC*TSTC - (TSTC*TSaxisV)*(TSTC*TSaxisV);
      if(TSdistTL[TSnlog]<0) TSdistTL[TSnlog] = 0;
      TSdistTL[TSnlog] = sqrt(TSdistTL[TSnlog]);
      //
      // computing TSdistT : using only the transverse position measurement
      //
      idents::CalXtalId id = recLog.getPackedId();
      itow = id.getTower();
      ilay = id.getLayer();
      icol = id.getColumn();
      itowy = itow/4;
      itowx = itow-4*itowy;
      if(ilay%2==0)
	{
	  TSxtalV = Vector(1,0,0);
	  TSxtalC = Point(-1.5*374.5+374.5*(double)itowx,TSxtalP.y(),TSxtalP.z());
	}
      else
	{
	  TSxtalV = Vector(0,1,0);
	  if(m_xNum==4 && m_yNum==1) // beamtest configuration
	    TSxtalC = Point(TSxtalP.x(),0,TSxtalP.z());
	  else
	    TSxtalC = Point(TSxtalP.x(),-1.5*374.5+374.5*(double)itowy,TSxtalP.z());
	}
      TSTC = TSxtalC-TSaxisP;
      lambda = 1-(TSxtalV*TSaxisV)*(TSxtalV*TSaxisV);
      if(lambda==0) // xtal axis and main axis are parallel
	{
	  TSdistT[TSnlog] = TSTC*TSTC - (TSTC*TSaxisV)*(TSTC*TSaxisV);
	  if(TSdistT[TSnlog]<0) TSdistT[TSnlog] = 0;
	  TSdistT[TSnlog] = sqrt(TSdistT[TSnlog]);
	}
      else
	{
	  lambda = (-(TSTC*TSxtalV)+(TSTC*TSaxisV)*(TSxtalV*TSaxisV))/lambda;
	  if(lambda>lambdamax)
	    lambda = lambdamax;
	  if(lambda<-lambdamax)
	    lambda = -lambdamax;
	  TSxtalP = TSxtalC + lambda*TSxtalV;
	  TSTC = TSxtalP-TSaxisP;
	  TSdistT[TSnlog] = TSTC*TSTC - (TSTC*TSaxisV)*(TSTC*TSaxisV);
	  if(TSdistT[TSnlog]<0) TSdistT[TSnlog] = 0;
	  TSdistT[TSnlog] = sqrt(TSdistT[TSnlog]);
	}
      ++TSnlog;
    }
  
  return 0;
}

int CalValsTool::TSfillTS(int optts)
{
  // optts = 0 : use only transverse position measurement
  // optts = 1 : use both transverse and longitudinal position measurements

  if(TSnlog<=0) return 1;
  int i,j;
  double TStotalenergy = 0;
  for(i=0;i<TSnlog;++i)
    {
      TSTS[i] = 0;
      TSiused[i] = 0;
      TSiorder[i] = -1;
      TStotalenergy += TSenergy[i];
      if(optts==0)
	TSdist[i] = TSdistT[i]; 
      else
	TSdist[i] = TSdistTL[i]; 
    }
  if(TStotalenergy<=0) return 1;
  double mindist = 9999999;
  int imin = 0;
  for(i=0;i<TSnlog;++i)
    {
      mindist = 9999999;
      for(j=0;j<TSnlog;++j)
	{
	  if(TSiused[j]) continue;
	  if(TSdist[j]<mindist)
	    {
	      mindist = TSdist[j];
	      imin = j;
	    }
	}
      TSiorder[i] = imin;
      TSiused[imin] = 1;
    }
  double efrac = 0;
  double weight = 0;
  double transversesize = 0;
  for(i=0;i<TSnlog;++i)
    {
      efrac += TSenergy[TSiorder[i]];
      weight += TSenergy[TSiorder[i]];
      transversesize += TSenergy[TSiorder[i]]*TSdist[TSiorder[i]]*TSdist[TSiorder[i]];
      TSefrac[i] = efrac/TStotalenergy;
      if(weight>0)
	TSTS[i] = sqrt(transversesize/weight);
    }

  return 0;
}

double CalValsTool::TSgetinterpolationTS(double efrac)
{
  if(TSnlog<=0) return -999;
  //
  if(efrac<=TSefrac[0]) return TSTS[0];
  if(efrac>=TSefrac[TSnlog-1]) return TSTS[TSnlog-1];
  //
  int i;
  int j = -1;
  for(i=0;i<TSnlog-1;++i)
    {
      if(TSefrac[i]<efrac && TSefrac[i+1]>=efrac)
	{
	  j = i;
	  break;
	}
    }
  if(j==-1) return TSTS[TSnlog-1];

  i = j;
  if(TSefrac[i+1]-TSefrac[i]>0)
    return ((TSefrac[i+1]-efrac)*TSTS[i]+(efrac-TSefrac[i])*TSTS[i+1])/(TSefrac[i+1]-TSefrac[i]);

  return (TSTS[i]+TSTS[i+1])/2;
}

void CalValsTool::zeroVals()
{
    ValBase::zeroVals();

    // This is so zeroing the Cal vals will keep CalEnergyRaw
    SmartDataPtr<Event::CalClusterCol>     
        pCals(m_pEventSvc,EventModel::CalRecon::CalClusterCol);
    if(pCals) {
        if (pCals->empty()) return;
        Event::CalCluster* calCluster = pCals->front();
        CAL_EnergyRaw  = calCluster->getCalParams().getEnergy();
    }
}
