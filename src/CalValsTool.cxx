/** @file CalValsTool.cxx
@brief Calculates the Cal analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/CalValsTool.cxx,v 1.138 2013/03/25 14:26:33 bruel Exp $
*/
//#define PRE_CALMOD 1

// Include files


#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/DataSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalParams.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalClusterMap.h"

#include "GaudiKernel/IToolSvc.h"
//#include "GlastSvc/Reco/IPropagatorTool.h"
//#include "GlastSvc/Reco/IPropagator.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrFlagHitsTool.h"

#include "CalUtil/IUBinterpolateTool.h"
#include "CalUtil/ICalClusterHitTool.h"

#include "idents/TowerId.h" 
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/Rotation.h"

#include "TMath.h"

#include "Doca.h"

/*! @class CaValsTool
@brief calculates Cal values

@authors Bill Atwood, Leon Rochester
*/

namespace {
    double truncFraction = 0.9;
    const int _badInt = -1;
    const float _badFloat = -2.0;
    double maxPath         = 2500.;  // limit the upward propagator
    double _vetoRegion     = 10.0;   // regions for Cal1SSDVeto hit counts

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

    /// pointer to queryclusterstool
    ITkrQueryClustersTool* pQueryClusters;
    /// pointer to flagHitsTool
    ITkrFlagHitsTool* pFlagHits;

  ICalClusterHitTool*       m_pCalClusterHitTool;

    /// some Geometry
    double m_towerPitch;
    int    m_xNum;
    int    m_yNum;
    int    m_nLayers;
    int    m_nCsI;

    /// gets the CAL info from detModel
    StatusCode getCalInfo();
    double activeDist(Point pos, int &view) const;
  
    /// SSD Veto function, emulated from Tkr one
    float CalSSDEvaluation(const Event::CalCluster* cluster);
    // Local variables to transfer results of SSD calculation
    float m_VetoPlaneCrossed; 
    float m_VetoTrials;
    float m_SSDVeto;
    float m_SSDVetoNoHitFlag;
    float m_VetoUnknown;
    float m_VetoDeadPlane;
    float m_VetoTruncated; 
    float m_VetoTower;   
    float m_VetoGapCorner;
    float m_VetoGapEdge;   
    float m_VetoBadCluster;
    float m_VetoHitFound;
    // flag for additional veto info
    bool  m_enableVetoDiagnostics;
    // properties
    double m_minVetoError;
    double m_maxVetoError;
    double m_vetoNSigma;
    double m_vetoDTheta;

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
    float CAL_TwrEdgeCntrCor;
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

    float CAL_MIP_Diff; 
    float CAL_MIP_Ratio;

    
    // Stuff to study clustering
    float CAL_NumXtals_uber;
    float CAL_NumTruncXtals_uber;
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

    float CAL_MomTransRms_uber;
    float CAL_MomLongRms_uber;
    float CAL_ClassGamProb_uber;
    float CAL_ClassHadProb_uber;
    float CAL_ClassGhostProb_uber;
    float CAL_ClassMipProb_uber;

    float CAL_num_clusters;
    float CAL_rest_energy;
    float CAL_rest_numXtals;

    // Variables referring to the first cluster---added after the restructuring
    // of the CAL moments analysis (Luca Baldini, Dec. 26, 2010).
    // Basic variables.
    float CAL_Clu1_NumXtals;
    float CAL_Clu1_NumTruncXtals;
    float CAL_Clu1_NumSaturatedXtals;
    float CAL_Clu1_RawEnergySum;
    float CAL_Clu1_CorrEnergySum;
    float CAL_Clu1_XtalEneMax;
    float CAL_Clu1_XtalEneRms;
    float CAL_Clu1_XtalEneSkewness;
    // Variables from the moments analysis.
    float CAL_Clu1_MomXCntr;
    float CAL_Clu1_MomYCntr;
    float CAL_Clu1_MomZCntr;
    float CAL_Clu1_MomXCntrCor;
    float CAL_Clu1_MomYCntrCor;
    float CAL_Clu1_MomZCntrCor;
    float CAL_Clu1_MomXDir;
    float CAL_Clu1_MomYDir;
    float CAL_Clu1_MomZDir;
    float CAL_Clu1_MomNumIterations;
  float CAL_Clu1_MomRunFlag; // FOR DEBUG ONLY!!!!
    float CAL_Clu1_MomNumCoreXtals;
    float CAL_Clu1_MomTransRms;
    float CAL_Clu1_MomLongRms;
    float CAL_Clu1_MomLongRmsAysm;
    float CAL_Clu1_MomLongSkewness;
    float CAL_Clu1_MomCoreEneFrac;
    float CAL_Clu1_MomFullLenght;
    float CAL_Clu1_MomdEdxAve;
    float CAL_Clu1_Mom_Cov_Cxx; // centroid
    float CAL_Clu1_Mom_Cov_Cyy;
    float CAL_Clu1_Mom_Cov_Czz;
    float CAL_Clu1_Mom_Cov_Dxx; // direction
    float CAL_Clu1_Mom_Cov_Dyy;
    float CAL_Clu1_Mom_Cov_Dzz;
    float CAL_Clu1_Mom_Cov_Dxy; 
    float CAL_Clu1_Mom_Cov_Dxz;
    float CAL_Clu1_Mom_Cov_Dyz;
    float CAL_Clu1_Mom_Cov_Sxx; // slope
    float CAL_Clu1_Mom_Cov_Syy;
    float CAL_Clu1_Mom_Cov_Sxy;  
  // Variables from Philippe's fit.
    float CAL_Clu1_FitXCntr;
    float CAL_Clu1_FitYCntr;
    float CAL_Clu1_FitZCntr;
    float CAL_Clu1_FitXDir;
    float CAL_Clu1_FitYDir;
    float CAL_Clu1_FitZDir;
    float CAL_Clu1_FitNumLayers;
    float CAL_Clu1_FitChiSquare;
    // Variables from the Minimum Spanning Tree clustering.
    float CAL_Clu1_MstMinEdgeLen;
    float CAL_Clu1_MstMaxEdgeLen;
    float CAL_Clu1_MstAveEdgeLen;
    float CAL_Clu1_MstRmsEdgeLen;
    float CAL_Clu1_MstAveTruncEdgeLen;
    float CAL_Clu1_MstRmsTruncEdgeLen;
    // Variables from the cluster classification.
    float CAL_Clu1_ClassGamProb;
    float CAL_Clu1_ClassHadProb;
    float CAL_Clu1_ClassGhostProb;
    float CAL_Clu1_ClassMipProb;
  // Variables fron SSD Veto
  float CAL_Clu1_SSDVeto;
  float CAL_Clu1_SSDVetoNoHitFlag;
  float CAL_Clu1_VetoPlaneCrossed;
  float CAL_Clu1_VetoHitFound;
  float CAL_Clu1_VetoTrials;
  float CAL_Clu1_VetoUnknown;
  float CAL_Clu1_VetoDeadPlane;
  float CAL_Clu1_VetoTruncated;
  float CAL_Clu1_VetoTower;     
  float CAL_Clu1_VetoGapCorner;
  float CAL_Clu1_VetoGapEdge;
  float CAL_Clu1_VetoBadCluster; 

  float CAL_Clu2_NumXtals;
  float CAL_Clu2_RawEnergySum;
  float CAL_Clu2_XtalEneMax;
  float CAL_Clu2_MomXCntr;
  float CAL_Clu2_MomYCntr;
  float CAL_Clu2_MomZCntr;
  float CAL_Clu2_MomXDir;
  float CAL_Clu2_MomYDir;
  float CAL_Clu2_MomZDir;
  float CAL_Clu2_MomNumCoreXtals;
  float CAL_Clu2_MomTransRms;
  float CAL_Clu2_MomLongRms;
  float CAL_Clu2_ClassGamProb;
  float CAL_Clu2_ClassHadProb;
  float CAL_Clu2_ClassGhostProb;
  float CAL_Clu2_ClassMipProb;
  float CAL_Clu2_Dist;
  float CAL_Clu2_Doca;

  float CAL_Clu3_NumXtals;
  float CAL_Clu3_RawEnergySum;
  float CAL_Clu3_XtalEneMax;
  float CAL_Clu3_MomXCntr;
  float CAL_Clu3_MomYCntr;
  float CAL_Clu3_MomZCntr;
  float CAL_Clu3_MomXDir;
  float CAL_Clu3_MomYDir;
  float CAL_Clu3_MomZDir;
  float CAL_Clu3_MomNumCoreXtals;
  float CAL_Clu3_MomTransRms;
  float CAL_Clu3_MomLongRms;
  float CAL_Clu3_ClassGamProb;
  float CAL_Clu3_ClassHadProb;
  float CAL_Clu3_ClassGhostProb;
  float CAL_Clu3_ClassMipProb;
  float CAL_Clu3_Dist;
  float CAL_Clu3_Doca;

  float CAL_Clu4_NumXtals;
  float CAL_Clu4_RawEnergySum;
  float CAL_Clu4_XtalEneMax;
  float CAL_Clu4_MomXCntr;
  float CAL_Clu4_MomYCntr;
  float CAL_Clu4_MomZCntr;
  float CAL_Clu4_MomXDir;
  float CAL_Clu4_MomYDir;
  float CAL_Clu4_MomZDir;
  float CAL_Clu4_MomNumCoreXtals;
  float CAL_Clu4_MomTransRms;
  float CAL_Clu4_MomLongRms;
  float CAL_Clu4_ClassGamProb;
  float CAL_Clu4_ClassHadProb;
  float CAL_Clu4_ClassGhostProb;
  float CAL_Clu4_ClassMipProb;
  float CAL_Clu4_Dist;
  float CAL_Clu4_Doca;

  float CAL_Clu5_NumXtals;
  float CAL_Clu5_RawEnergySum;
  float CAL_Clu5_XtalEneMax;
  float CAL_Clu5_MomXCntr;
  float CAL_Clu5_MomYCntr;
  float CAL_Clu5_MomZCntr;
  float CAL_Clu5_MomXDir;
  float CAL_Clu5_MomYDir;
  float CAL_Clu5_MomZDir;
  float CAL_Clu5_MomNumCoreXtals;
  float CAL_Clu5_MomTransRms;
  float CAL_Clu5_MomLongRms;
  float CAL_Clu5_ClassGamProb;
  float CAL_Clu5_ClassHadProb;
  float CAL_Clu5_ClassGhostProb;
  float CAL_Clu5_ClassMipProb;
  float CAL_Clu5_Dist;
  float CAL_Clu5_Doca;

  float CAL_Uber2_NumXtals;
  float CAL_Uber2_RawEnergySum;
  float CAL_Uber2_XtalEneMax;
  float CAL_Uber2_MomXCntr;
  float CAL_Uber2_MomYCntr;
  float CAL_Uber2_MomZCntr;
  float CAL_Uber2_MomXDir;
  float CAL_Uber2_MomYDir;
  float CAL_Uber2_MomZDir;
  float CAL_Uber2_MomNumCoreXtals;
  float CAL_Uber2_MomTransRms;
  float CAL_Uber2_MomLongRms;



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

  float CAL_newcfp_energy;      // Energy from Full Profile tool
  float CAL_newcfp_energy_err;
  float CAL_newcfp_totChiSq;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_seltotChiSq;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_ChiSqdist;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_Chidist;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_energyUB;    // Energy from Full Profile tool unbiased
  float CAL_newcfp_calEffRLn;   // Effective radiation lengths in the Cal
  float CAL_newcfp_tkrRLn;   // Effective radiation lengths in the tkr
  float CAL_newcfp_par0;
  float CAL_newcfp_par1;
  float CAL_newcfp_parcf;
  float CAL_newcfp_xcor;
  float CAL_newcfp_ycor;
  float CAL_newcfp_alpha;       // fit parameter alpha
  float CAL_newcfp_tmax;        // fit parameter tmax
  float CAL_newcfp_fiterrflg;   // fit error flag
  float CAL_newcfp_widening;   // 
  float CAL_newcfp_nxtalsat;
  float CAL_newcfp_nxtalsel;
  float CAL_newcfp_nlayersat;

  float CAL_newcfp_p0;
  float CAL_newcfp_p1;
  float CAL_newcfp_p2;
  float CAL_newcfp_v0;
  float CAL_newcfp_v1;
  float CAL_newcfp_v2;

  float CAL_newcfp_calfit_energy;      // Energy from Full Profile tool
  float CAL_newcfp_calfit_energyub;    // Energy ub version for cal-only
  float CAL_newcfp_calfit_energy_err;
  float CAL_newcfp_calfit_totChiSq;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_calfit_seltotChiSq;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_calfit_ChiSqdist;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_calfit_Chidist;    // Total ChiSquare of fit divided by 11
  float CAL_newcfp_calfit_calEffRLn;   // Effective radiation lengths in the Cal
  float CAL_newcfp_calfit_par0;
  float CAL_newcfp_calfit_par1;
  float CAL_newcfp_calfit_parcf;
  float CAL_newcfp_calfit_xcor;
  float CAL_newcfp_calfit_ycor;
  float CAL_newcfp_calfit_alpha;       // fit parameter alpha
  float CAL_newcfp_calfit_tmax;        // fit parameter tmax
  float CAL_newcfp_calfit_fiterrflg;   // fit error flag
  float CAL_newcfp_calfit_widening;   // 
  float CAL_newcfp_calfit_nxtalsel;

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
    float CAL_Track_DOCA_Cor;
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

    float CAL_uber_TS_TKR_T_100;
    float CAL_uber_TS_TKR_TL_100;
    float CAL_uber2_TS_TKR_T_100;
    float CAL_uber2_TS_TKR_TL_100;

    int TSfillTSdist(std::vector<Event::CalXtalRecData*> xtallist);
    int TSfillTS(int optts);
    double TSgetinterpolationTS(double efrac);

    IUBinterpolateTool* m_ubInterpolateTool;

};

namespace {
    // this is the test distance for the CAL_EdgeEnergy variable
    double _deltaEdge = 50.0;
}

// Static factory for instantiation of algtool objects
//static ToolFactory<CalValsTool> s_factory;
//const IToolFactory& CalValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(CalValsTool);

// Standard Constructor
CalValsTool::CalValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 

    // params for SSD Veto
    declareProperty("enableVetoDiagnostics", m_enableVetoDiagnostics=true);
    declareProperty("minVetoError", m_minVetoError=1.);
    declareProperty("maxVetoError", m_maxVetoError=100000.0);
    declareProperty("vetoNSigma", m_vetoNSigma=2.0);
    declareProperty("vetoDTheta", m_vetoDTheta=0.03); // ONLY for simplified error matrix


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
<tr><td> CalTwrEdgeCntrCor
<td>F<td>   Distance of the energy centroid (corrected for the hodoscopic effect) from the nearest tower boundary.  
<tr><td> CalGapFraction 
<td>F<td>   Approximate fraction of the shower volumn which falls in inter-tower gaps. 
<tr><td> CalTrackSep 
<td>F<td>   Distance between impact points of two best tracks at CAL front face; 
zero if only one track 
<tr><td> CalTrackDoca 
<td>F<td>   Distance between the projected vertex (or track if only one track) 
and the energy centroid, evaluated at the z of the centroid. 
<tr><td> CalTrackAngle 
<td>F<td>   Angle between TKR track direction and CAL cluster direction
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
<tr><td> CalCfpEnergyUB
<td>F<td> unbiased CalCfpEnergy
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
<tr><td> Cal11SSDVeto  
<td>F<td> A version of SSD Veto that uses Cal moments analysis direction. Emulating Tkr calculation.
<tr><td> Cal1SSDVetoPlaneCrossed
<td>I<td>   Number of planes contributiong to the SSD Veto. This doesn't count the points
where a track crosses in a gap.
<tr><td> Cal1MomCov[CXX,CYY,CZZ, DXX,DYY,DZZ, DXY,DXZ,DYZ, SXX,SYY,SXY ] 
<td>F<td> Non-zero element of Cal covariance matrix : C = is for centroid elements, D for direction (6 elements), S for slope (3 elements)

<tr><td> CalNewCfpEnergy
<td>F<td> DOCUMENTATION NEEDED
<tr><td> CalNewCfpSelChiSq 
<td>F<td> DOCUMENTATION NEEDED

<tr><td> Cal1NumXtals
<td>F<td> Number of xtals in the first CAL cluster

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

        if(!toolSvc->retrieveTool("TkrQueryClustersTool", pQueryClusters)) {
          log << MSG::ERROR << "Couldn't retrieve TkrQueryClusterTool" << endreq;
          return StatusCode::FAILURE;
        }
        
        if(!toolSvc->retrieveTool("TkrFlagHitsTool", pFlagHits)) {
          log << MSG::ERROR << "Couldn't retrieve TkrFlagHitsTool" << endreq;
          return StatusCode::FAILURE;
        }

        if(!toolSvc->retrieveTool("UBinterpolateTool", m_ubInterpolateTool)) {
          log << MSG::ERROR << "Couldn't retrieve UBinterpolateTool" << endreq;
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

        m_pCalClusterHitTool = 0;
        sc = pToolSvc->retrieveTool("CalClusterHitTool", m_pCalClusterHitTool);
        if( sc.isFailure() ) {
          log << MSG::ERROR << "Unable to find tool: " "CalClusterHitTool" << endreq;
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

    addItem("CalEnergyRaw",  &CAL_EnergyRaw, true);
    addItem("CalEnergyCorr", &CAL_EnergyCorr);

    addItem("CalLeakCorr",   &CAL_Leak_Corr); 
    addItem("CalEdgeCorr",   &CAL_Edge_Corr, true);
    addItem("CalTotalCorr",  &CAL_Total_Corr);

    addItem("CalCsIRLn",     &CAL_CsI_RLn, true);
    addItem("CalTotRLn",     &CAL_Tot_RLn);
    addItem("CalCntRLn",     &CAL_Cntr_RLn);
    addItem("CalLATRLn",     &CAL_LAT_RLn, true);
    addItem("CalDeadTotRat", &CAL_DeadTot_Rat);
    addItem("CalDeadCntRat", &CAL_DeadCnt_Rat);

    addItem("CalTPred",      &CAL_t_Pred);
    addItem("CalDeltaT",     &CAL_deltaT);

    addItem("CalGapFraction",&CAL_Gap_Fraction, true);
    addItem("CalTwrEdgeCntr",&CAL_TwrEdgeCntr, true);
    addItem("CalTwrEdgeCntrCor",&CAL_TwrEdgeCntrCor);
    addItem("CalTwrEdge",    &CAL_TwrEdgeTop, true);
    addItem("CalLATEdge",    &CAL_LATEdge,    true);
    addItem("CalEdgeEnergy", &CAL_EdgeEnergy, true);
    addItem("CalTrackDoca",  &CAL_Track_DOCA, true);
    addItem("CalTrackDocacor",  &CAL_Track_DOCA_Cor, true);
    addItem("CalTrackDocaCor",  &CAL_Track_DOCA_Cor);
    addItem("CalTrackAngle", &CAL_Track_Angle, true);
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
    addItem("CalBkHalfRatio",&CAL_BkHalf_Ratio, true);

    addItem("CalXtalsTrunc", &CAL_Num_Xtals_Trunc);
    addItem("CalNumXtals",   &CAL_Num_Xtals);
    addItem("CalXtalRatio",  &CAL_Xtal_Ratio);
    addItem("CalXtalMaxEne", &CAL_Xtal_maxEne, true);
    addItem("CalMaxNumXtalsInLayer", 
        &CAL_Max_Num_Xtals_In_Layer);

    addItem("CalTransRms",    &CAL_Trans_Rms, true);
    addItem("CalLongRms",     &CAL_Long_Rms,  true);
    addItem("CalLRmsAsym",    &CAL_LRms_Asym);

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
    addItem("CalTopGapDist", &CAL_Top_Gap_Dist, true);

    addItem("CalXEcntr2",     &CAL_xEcntr2);
    addItem("CalYEcntr2",     &CAL_yEcntr2);
    addItem("CalZEcntr2",     &CAL_zEcntr2);
    addItem("CalXDir2",       &CAL_xdir2);
    addItem("CalYDir2",       &CAL_ydir2);
    addItem("CalZDir2",       &CAL_zdir2);
    addItem("CalPosDirChisq",       &CAL_posdir_chisq);
    addItem("CalPosDirNLayers",       &CAL_posdir_nlayers);
    addItem("CalNSaturated",       &CAL_nsaturated, true);

    addItem("CalUberNumXtals",  &CAL_NumXtals_uber);
    addItem("CalUberNumTruncXtals",  &CAL_NumTruncXtals_uber);
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

    addItem("CalUberTransRms",    &CAL_MomTransRms_uber);
    addItem("CalUberLongRms",     &CAL_MomLongRms_uber);
    addItem("CalUberGamProb",     &CAL_ClassGamProb_uber);
    addItem("CalUberHadProb",     &CAL_ClassHadProb_uber);
    addItem("CalUberGhostProb",   &CAL_ClassGhostProb_uber);
    addItem("CalUberMipProb",     &CAL_ClassMipProb_uber);

    addItem("CalNumClusters",     &CAL_num_clusters);
    addItem("CalRestEnergy",      &CAL_rest_energy);
    addItem("CalRestNumXtals",    &CAL_rest_numXtals);

    addItem("CalTrkXtalRms",       &CAL_track_rms, true);
    addItem("CalTrkXtalRmsE",      &CAL_track_E_rms, true);
    addItem("CalTrkXtalRmsTrunc",  &CAL_track_rms_trunc);
    addItem("CalTrkXtalRmsETrunc", &CAL_track_E_rms_trunc, true);

    // Variables referring to the first cluster---added after the restructuring
    // of the CAL moments analysis (Luca Baldini, Dec. 26, 2010).
    // Basic variables.
    addItem("Cal1NumXtals",  &CAL_Clu1_NumXtals, true);
    addItem("Cal1NumTruncXtals",  &CAL_Clu1_NumTruncXtals, true);
    addItem("Cal1NumSaturatedXtals",  &CAL_Clu1_NumSaturatedXtals, true);
    addItem("Cal1RawEnergySum",  &CAL_Clu1_RawEnergySum);
    addItem("Cal1CorrEnergySum",  &CAL_Clu1_CorrEnergySum);
    addItem("Cal1XtalEneMax",  &CAL_Clu1_XtalEneMax, true);
    addItem("Cal1XtalEneRms",  &CAL_Clu1_XtalEneRms, true);
    addItem("Cal1XtalEneSkewness",  &CAL_Clu1_XtalEneSkewness);
    // Variables from the moments analysis.
    addItem("Cal1MomXCntr",  &CAL_Clu1_MomXCntr);
    addItem("Cal1MomYCntr",  &CAL_Clu1_MomYCntr);
    addItem("Cal1MomZCntr",  &CAL_Clu1_MomZCntr);
    addItem("Cal1MomXCntrcor",  &CAL_Clu1_MomXCntrCor);
    addItem("Cal1MomYCntrcor",  &CAL_Clu1_MomYCntrCor);
    addItem("Cal1MomZCntrcor",  &CAL_Clu1_MomZCntrCor);
    addItem("Cal1MomXCntrCor",  &CAL_Clu1_MomXCntrCor);
    addItem("Cal1MomYCntrCor",  &CAL_Clu1_MomYCntrCor);
    addItem("Cal1MomZCntrCor",  &CAL_Clu1_MomZCntrCor);
    addItem("Cal1MomXDir",  &CAL_Clu1_MomXDir, true);
    addItem("Cal1MomYDir",  &CAL_Clu1_MomYDir, true);
    addItem("Cal1MomZDir",  &CAL_Clu1_MomZDir, true);
    addItem("Cal1MomNumIterations",  &CAL_Clu1_MomNumIterations);
    addItem("Cal1MomRunFlag",  &CAL_Clu1_MomRunFlag); // FOR DEBUG ONLY!!!!
    addItem("Cal1MomNumCoreXtals",  &CAL_Clu1_MomNumCoreXtals, true);
    addItem("Cal1TransRms",  &CAL_Clu1_MomTransRms, true);
    addItem("Cal1LongRms",  &CAL_Clu1_MomLongRms, true);
    addItem("Cal1LongRmsAsym",  &CAL_Clu1_MomLongRmsAysm);
    addItem("Cal1LongSkewness",  &CAL_Clu1_MomLongSkewness);
    addItem("Cal1CoreEneFrac",  &CAL_Clu1_MomCoreEneFrac);
    addItem("Cal1FullLength",  &CAL_Clu1_MomFullLenght);
    addItem("Cal1dEdxAve",  &CAL_Clu1_MomdEdxAve);
    // From Cov Matrix - cntr + slope
    addItem("Cal1MomCovCXX",  &CAL_Clu1_Mom_Cov_Cxx);// centroid
    addItem("Cal1MomCovCYY",  &CAL_Clu1_Mom_Cov_Cyy);
    addItem("Cal1MomCovCZZ",  &CAL_Clu1_Mom_Cov_Czz);
    addItem("Cal1MomCovDXX",  &CAL_Clu1_Mom_Cov_Dxx);// direction
    addItem("Cal1MomCovDYY",  &CAL_Clu1_Mom_Cov_Dyy);
    addItem("Cal1MomCovDZZ",  &CAL_Clu1_Mom_Cov_Dzz);
    addItem("Cal1MomCovDXY",  &CAL_Clu1_Mom_Cov_Dxy);
    addItem("Cal1MomCovDXZ",  &CAL_Clu1_Mom_Cov_Dxz);
    addItem("Cal1MomCovDYZ",  &CAL_Clu1_Mom_Cov_Dyz);
    addItem("Cal1MomCovSXX",  &CAL_Clu1_Mom_Cov_Sxx);// slope
    addItem("Cal1MomCovSYY",  &CAL_Clu1_Mom_Cov_Syy);
    addItem("Cal1MomCovSXY",  &CAL_Clu1_Mom_Cov_Sxy);
    // Variables from Philippe's fit.
    addItem("Cal1FitXCntr",  &CAL_Clu1_FitXCntr);
    addItem("Cal1FitYCntr",  &CAL_Clu1_FitYCntr);
    addItem("Cal1FitZCntr",  &CAL_Clu1_FitZCntr);
    addItem("Cal1FitXDir",  &CAL_Clu1_FitXDir);
    addItem("Cal1FitYDir",  &CAL_Clu1_FitYDir);
    addItem("Cal1FitZDir",  &CAL_Clu1_FitZDir);
    addItem("Cal1FitNumLayers",  &CAL_Clu1_FitNumLayers);
    addItem("Cal1FitChiSquare",  &CAL_Clu1_FitChiSquare, true);
    // Variables from the Minimum Spanning Tree clustering.
    addItem("Cal1MstMinEdgeLen",  &CAL_Clu1_MstMinEdgeLen);
    addItem("Cal1MstMaxEdgeLen",  &CAL_Clu1_MstMaxEdgeLen);
    addItem("Cal1MstAveEdgeLen",  &CAL_Clu1_MstAveEdgeLen);
    addItem("Cal1MstRmsEdgeLen",  &CAL_Clu1_MstRmsEdgeLen);
    addItem("Cal1MstAveTruncEdgeLen",  &CAL_Clu1_MstAveTruncEdgeLen);
    addItem("Cal1MstRmsTruncEdgeLen",  &CAL_Clu1_MstRmsTruncEdgeLen);
    // Variables from the cluster classification.
    addItem("Cal1GamProb",  &CAL_Clu1_ClassGamProb);
    addItem("Cal1HadProb",  &CAL_Clu1_ClassHadProb);
    addItem("Cal1GhostProb",  &CAL_Clu1_ClassGhostProb);
    addItem("Cal1MipProb",  &CAL_Clu1_ClassMipProb);
    // Variables for SSD Veto
    addItem("Cal1SSDVeto",  &CAL_Clu1_SSDVeto);
    addItem("Cal1SSDVetoNoHitFlag",  &CAL_Clu1_SSDVetoNoHitFlag);
    addItem("Cal1SSDVetoPlaneCrossed",  &CAL_Clu1_VetoPlaneCrossed);

    if(m_enableVetoDiagnostics) {
        addItem("Cal1SSDVetoHitFound",   &CAL_Clu1_VetoHitFound);
        addItem("Cal1SSDVetoTrials",     &CAL_Clu1_VetoTrials);
        addItem("Cal1SSDVetoUnknown",    &CAL_Clu1_VetoUnknown);
        addItem("Cal1SSDVetoTower",      &CAL_Clu1_VetoTower);
        addItem("Cal1SSDVetoGapCorner",  &CAL_Clu1_VetoGapCorner);
        addItem("Cal1SSDVetoGapEdge",    &CAL_Clu1_VetoGapEdge);
        addItem("Cal1SSDVetoBadCluster", &CAL_Clu1_VetoBadCluster);
        addItem("Cal1SSDVetoDeadPlane",  &CAL_Clu1_VetoDeadPlane);
        addItem("Cal1SSDVetoTruncated",  &CAL_Clu1_VetoTruncated);
    }

    addItem("Cal2NumXtals",  &CAL_Clu2_NumXtals);
    addItem("Cal2RawEnergySum",  &CAL_Clu2_RawEnergySum);
    addItem("Cal2XtalEneMax",  &CAL_Clu2_XtalEneMax);
    addItem("Cal2MomXCntr",  &CAL_Clu2_MomXCntr);
    addItem("Cal2MomYCntr",  &CAL_Clu2_MomYCntr);
    addItem("Cal2MomZCntr",  &CAL_Clu2_MomZCntr);
    addItem("Cal2MomXDir",  &CAL_Clu2_MomXDir);
    addItem("Cal2MomYDir",  &CAL_Clu2_MomYDir);
    addItem("Cal2MomZDir",  &CAL_Clu2_MomZDir);
    addItem("Cal2MomNumCoreXtals",  &CAL_Clu2_MomNumCoreXtals);
    addItem("Cal2TransRms",  &CAL_Clu2_MomTransRms);
    addItem("Cal2LongRms",  &CAL_Clu2_MomLongRms);
    addItem("Cal2GamProb",  &CAL_Clu2_ClassGamProb);
    addItem("Cal2HadProb",  &CAL_Clu2_ClassHadProb);
    addItem("Cal2GhostProb",  &CAL_Clu2_ClassGhostProb);
    addItem("Cal2MipProb",  &CAL_Clu2_ClassMipProb);
    addItem("Cal2Dist",  &CAL_Clu2_Dist);
    addItem("Cal2Doca",  &CAL_Clu2_Doca);

    addItem("Cal3NumXtals",  &CAL_Clu3_NumXtals);
    addItem("Cal3RawEnergySum",  &CAL_Clu3_RawEnergySum);
    addItem("Cal3XtalEneMax",  &CAL_Clu3_XtalEneMax);
    addItem("Cal3MomXCntr",  &CAL_Clu3_MomXCntr);
    addItem("Cal3MomYCntr",  &CAL_Clu3_MomYCntr);
    addItem("Cal3MomZCntr",  &CAL_Clu3_MomZCntr);
    addItem("Cal3MomXDir",  &CAL_Clu3_MomXDir);
    addItem("Cal3MomYDir",  &CAL_Clu3_MomYDir);
    addItem("Cal3MomZDir",  &CAL_Clu3_MomZDir);
    addItem("Cal3MomNumCoreXtals",  &CAL_Clu3_MomNumCoreXtals);
    addItem("Cal3TransRms",  &CAL_Clu3_MomTransRms);
    addItem("Cal3LongRms",  &CAL_Clu3_MomLongRms);
    addItem("Cal3GamProb",  &CAL_Clu3_ClassGamProb);
    addItem("Cal3HadProb",  &CAL_Clu3_ClassHadProb);
    addItem("Cal3GhostProb",  &CAL_Clu3_ClassGhostProb);
    addItem("Cal3MipProb",  &CAL_Clu3_ClassMipProb);
    addItem("Cal3Dist",  &CAL_Clu3_Dist);
    addItem("Cal3Doca",  &CAL_Clu3_Doca);

    addItem("Cal4NumXtals",  &CAL_Clu4_NumXtals);
    addItem("Cal4RawEnergySum",  &CAL_Clu4_RawEnergySum);
    addItem("Cal4XtalEneMax",  &CAL_Clu4_XtalEneMax);
    addItem("Cal4MomXCntr",  &CAL_Clu4_MomXCntr);
    addItem("Cal4MomYCntr",  &CAL_Clu4_MomYCntr);
    addItem("Cal4MomZCntr",  &CAL_Clu4_MomZCntr);
    addItem("Cal4MomXDir",  &CAL_Clu4_MomXDir);
    addItem("Cal4MomYDir",  &CAL_Clu4_MomYDir);
    addItem("Cal4MomZDir",  &CAL_Clu4_MomZDir);
    addItem("Cal4MomNumCoreXtals",  &CAL_Clu4_MomNumCoreXtals);
    addItem("Cal4TransRms",  &CAL_Clu4_MomTransRms);
    addItem("Cal4LongRms",  &CAL_Clu4_MomLongRms);
    addItem("Cal4GamProb",  &CAL_Clu4_ClassGamProb);
    addItem("Cal4HadProb",  &CAL_Clu4_ClassHadProb);
    addItem("Cal4GhostProb",  &CAL_Clu4_ClassGhostProb);
    addItem("Cal4MipProb",  &CAL_Clu4_ClassMipProb);
    addItem("Cal4Dist",  &CAL_Clu4_Dist);
    addItem("Cal4Doca",  &CAL_Clu4_Doca);

    addItem("Cal5NumXtals",  &CAL_Clu5_NumXtals);
    addItem("Cal5RawEnergySum",  &CAL_Clu5_RawEnergySum);
    addItem("Cal5XtalEneMax",  &CAL_Clu5_XtalEneMax);
    addItem("Cal5MomXCntr",  &CAL_Clu5_MomXCntr);
    addItem("Cal5MomYCntr",  &CAL_Clu5_MomYCntr);
    addItem("Cal5MomZCntr",  &CAL_Clu5_MomZCntr);
    addItem("Cal5MomXDir",  &CAL_Clu5_MomXDir);
    addItem("Cal5MomYDir",  &CAL_Clu5_MomYDir);
    addItem("Cal5MomZDir",  &CAL_Clu5_MomZDir);
    addItem("Cal5MomNumCoreXtals",  &CAL_Clu5_MomNumCoreXtals);
    addItem("Cal5TransRms",  &CAL_Clu5_MomTransRms);
    addItem("Cal5LongRms",  &CAL_Clu5_MomLongRms);
    addItem("Cal5GamProb",  &CAL_Clu5_ClassGamProb);
    addItem("Cal5HadProb",  &CAL_Clu5_ClassHadProb);
    addItem("Cal5GhostProb",  &CAL_Clu5_ClassGhostProb);
    addItem("Cal5MipProb",  &CAL_Clu5_ClassMipProb);
    addItem("Cal5Dist",  &CAL_Clu5_Dist);
    addItem("Cal5Doca",  &CAL_Clu5_Doca);

    addItem("CalUber2NumXtals",  &CAL_Uber2_NumXtals);
    addItem("CalUber2RawEnergySum",  &CAL_Uber2_RawEnergySum);
    addItem("CalUber2XtalEneMax",  &CAL_Uber2_XtalEneMax);
    addItem("CalUber2MomXCntr",  &CAL_Uber2_MomXCntr);
    addItem("CalUber2MomYCntr",  &CAL_Uber2_MomYCntr);
    addItem("CalUber2MomZCntr",  &CAL_Uber2_MomZCntr);
    addItem("CalUber2MomXDir",  &CAL_Uber2_MomXDir);
    addItem("CalUber2MomYDir",  &CAL_Uber2_MomYDir);
    addItem("CalUber2MomZDir",  &CAL_Uber2_MomZDir);
    addItem("CalUber2MomNumCoreXtals",  &CAL_Uber2_MomNumCoreXtals);
    addItem("CalUber2TransRms",  &CAL_Uber2_MomTransRms);
    addItem("CalUber2LongRms",  &CAL_Uber2_MomLongRms);

    addItem("CalCfpEnergy",  &CAL_cfp_energy, true);
    addItem("CalCfpChiSq",   &CAL_cfp_totChiSq);
    addItem("CalCfpEnergyUB",&CAL_cfp_energyUB, true);
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

    addItem("CalNewCfpEnergy",  &CAL_newcfp_energy, true);
    addItem("CalNewCfpEnergyErr",  &CAL_newcfp_energy_err);
    addItem("CalNewCfpChiSq",   &CAL_newcfp_totChiSq, true);
    addItem("CalNewCfpSelChiSq",   &CAL_newcfp_seltotChiSq, true);
    addItem("CalNewCfpSelChiSqDist",   &CAL_newcfp_ChiSqdist);
    addItem("CalNewCfpSelChiDist",   &CAL_newcfp_Chidist);
    addItem("CalNewCfpEnergyUB",&CAL_newcfp_energyUB, true);
    addItem("CalNewCfpEffRLn",  &CAL_newcfp_calEffRLn);
    addItem("CalNewCfpTkrRLn",  &CAL_newcfp_tkrRLn);
    addItem("CalNewCfpPar0",  &CAL_newcfp_par0);
    addItem("CalNewCfpPar1",  &CAL_newcfp_par1);
    addItem("CalNewCfpParcf",  &CAL_newcfp_parcf);
    addItem("CalNewCfpXcor",  &CAL_newcfp_xcor);
    addItem("CalNewCfpYcor",  &CAL_newcfp_ycor);
    addItem("CalNewCfpAlpha",  &CAL_newcfp_alpha);
    addItem("CalNewCfpTmax",  &CAL_newcfp_tmax);
    addItem("CalNewCfpFitErrFlg",  &CAL_newcfp_fiterrflg);
    addItem("CalNewCfpNLsat",  &CAL_newcfp_nlayersat);
    addItem("CalNewCfpNsat",  &CAL_newcfp_nxtalsat);
    addItem("CalNewCfpNsel",  &CAL_newcfp_nxtalsel);
    addItem("CalNewCfpPx",  &CAL_newcfp_p0);
    addItem("CalNewCfpPy",  &CAL_newcfp_p1);
    addItem("CalNewCfpPz",  &CAL_newcfp_p2);
    addItem("CalNewCfpVx",  &CAL_newcfp_v0, true);
    addItem("CalNewCfpVy",  &CAL_newcfp_v1, true);
    addItem("CalNewCfpVz",  &CAL_newcfp_v2, true);
    addItem("CalNewCfpCalEnergy",  &CAL_newcfp_calfit_energy, true);
    addItem("CalNewCfpCalEnergyUB",  &CAL_newcfp_calfit_energyub, true);
    addItem("CalNewCfpCalEnergyErr",  &CAL_newcfp_calfit_energy_err);
    addItem("CalNewCfpCalChiSq",   &CAL_newcfp_calfit_totChiSq, true);
    addItem("CalNewCfpCalSelChiSq",   &CAL_newcfp_calfit_seltotChiSq);
    addItem("CalNewCfpCalSelChiSqDist",   &CAL_newcfp_calfit_ChiSqdist);
    addItem("CalNewCfpCalSelChiDist",   &CAL_newcfp_calfit_Chidist);
    addItem("CalNewCfpCalEffRLn",  &CAL_newcfp_calfit_calEffRLn);
    addItem("CalNewCfpCalPar0",  &CAL_newcfp_calfit_par0);
    addItem("CalNewCfpCalPar1",  &CAL_newcfp_calfit_par1);
    addItem("CalNewCfpCalParcf",  &CAL_newcfp_calfit_parcf);
    addItem("CalNewCfpCalXcor",  &CAL_newcfp_calfit_xcor);
    addItem("CalNewCfpCalYcor",  &CAL_newcfp_calfit_ycor);
    addItem("CalNewCfpCalAlpha",  &CAL_newcfp_calfit_alpha);
    addItem("CalNewCfpCalTmax",  &CAL_newcfp_calfit_tmax, true);
    addItem("CalNewCfpCalFitErrFlg",  &CAL_newcfp_calfit_fiterrflg);
    addItem("CalNewCfpCalNsel",  &CAL_newcfp_calfit_nxtalsel);

    //addItem("CalLllEnergy",  &CAL_lll_energy);
    //addItem("CalLllEneErr",  &CAL_lll_energyErr);
    //addItem("CalTklEnergy",  &CAL_tkl_energy);
    //addItem("CalTklEneErr",  &CAL_tkl_energyErr);
    addItem("CalLkHdEnergy", &CAL_LkHd_energy);
    addItem("CalLkHdEneErr", &CAL_LkHd_energyErr);
    addItem("CalLkHdEnergyUB", &CAL_LkHd_energyUB);

    addItem("CalRmsLayerE",  &CAL_rmsLayerE, true);
    addItem("CalRmsLayerEBack",  &CAL_rmsLayerEBack, true);
    addItem("CalNLayersRmsBack", &CAL_nLayersRmsBack);
    addItem("CalEAveBack", &CAL_eAveBack);
    addItem("CalLayer0Ratio", &CAL_layer0Ratio, true);
    addItem("CalXPosRmsLL", &CAL_xPosRmsLastLayer);
    addItem("CalYPosRmsLL", &CAL_yPosRmsLastLayer);

    addItem("CalTrSizeCalT68",&CAL_TS_CAL_T_68);
    addItem("CalTrSizeCalT90",&CAL_TS_CAL_T_90);
    addItem("CalTrSizeCalT95",&CAL_TS_CAL_T_95, true);
    addItem("CalTrSizeCalT99",&CAL_TS_CAL_T_99);
    addItem("CalTrSizeCalT100",&CAL_TS_CAL_T_100, true);
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
    addItem("CalTrSizeTkrT100",&CAL_TS_TKR_T_100, true);
    addItem("CalTrSizeTkrTL68",&CAL_TS_TKR_TL_68);
    addItem("CalTrSizeTkrTL90",&CAL_TS_TKR_TL_90);
    addItem("CalTrSizeTkrTL95",&CAL_TS_TKR_TL_95);
    addItem("CalTrSizeTkrTL99",&CAL_TS_TKR_TL_99);
    addItem("CalTrSizeTkrTL100",&CAL_TS_TKR_TL_100);

    addItem("CalUberTrSizeTkrT100",&CAL_uber_TS_TKR_T_100);
    addItem("CalUberTrSizeTkrTL100",&CAL_uber_TS_TKR_TL_100);
    addItem("CalUber2TrSizeTkrT100",&CAL_uber2_TS_TKR_T_100);
    addItem("CalUber2TrSizeTkrTL100",&CAL_uber2_TS_TKR_TL_100);

    zeroVals();

    // adding bias maps
    //m_ubInterpolateTool->addBiasMap("Prof",           "$(ANALYSISNTUPLEDATAPATH)/BiasMapCalCfpEnergy.txt");
    //m_ubInterpolateTool->addBiasMap("NewProfalOnly", "$(ANALYSISNTUPLEDATAPATH)/BiasMapCalNewCfpCalEnergy_max20.txt");

    m_ubInterpolateTool->addBiasMap("Prof",           "$(CALUTILXMLPATH)/BiasMapCalCfpEnergy.txt"); 
    m_ubInterpolateTool->addBiasMap("NewProfCalOnly", "$(CALUTILXMLPATH)/BiasMapCalNewCfpCalEnergy_max20.txt");
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
    SmartDataPtr<Event::CalClusterMap>
      pCalClusterMap(m_pEventSvc,EventModel::CalRecon::CalClusterMap); 

    if(!pCalClusterMap) return sc;

    Event::CalCluster* uberCluster = (*pCalClusterMap).getUberCluster();
    Event::CalCluster* uber2Cluster = (*pCalClusterMap).getUber2Cluster();
    Event::CalClusterVec rawClusterVec = (*pCalClusterMap).get(EventModel::CalRecon::CalRawClusterVec);
    Event::CalCluster* firstCluster = rawClusterVec.front();

    CAL_num_clusters  = rawClusterVec.size();

    if(uberCluster==NULL) return sc;
    if(firstCluster==NULL) return sc;

    m_pCalClusterHitTool->fillRecDataVec(firstCluster);
    std::vector<Event::CalXtalRecData*> firstClusterXtalList = m_pCalClusterHitTool->getRecDataVec();

    // Recover pointer to CalEventEnergy info  
#ifdef PRE_CALMOD
    Event::CalEventEnergy* calEventEnergy = 
        SmartDataPtr<Event::CalEventEnergy>(m_pEventSvc, EventModel::CalRecon::CalEventEnergy);
#else
    Event::CalEventEnergy* calEventEnergy = 0;
    Event::CalEventEnergyMap* calEventEnergyMap = 
        SmartDataPtr<Event::CalEventEnergyMap>(m_pEventSvc,EventModel::CalRecon::CalEventEnergyMap) ;

    if (calEventEnergyMap && !calEventEnergyMap->empty())
    {
        Event::CalEventEnergyMap::iterator calEnergyItr = calEventEnergyMap->find(firstCluster);

        if (calEnergyItr != calEventEnergyMap->end()) calEventEnergy = calEnergyItr->second.front();
    }
#endif

    // Get relation between clusters and xtals
    Event::CalClusterHitTabList* xTal2ClusTabList = SmartDataPtr<Event::CalClusterHitTabList>(m_pEventSvc,EventModel::CalRecon::CalClusterHitTab);
    Event::CalClusterHitTab* xTal2ClusTab = 0;
    if (xTal2ClusTabList) xTal2ClusTab = new Event::CalClusterHitTab(xTal2ClusTabList);


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
                    //  std::cout << "CalValsTool CAL_cfp_energy " << CAL_cfp_energy << ' ' << CAL_cfp_energyUB << std::endl;
                }
                else {
                    float tkr1ZDir = -1;
                    if ( m_pTkrTool->getValCheck("Tkr1ZDir", tkr1ZDir).isSuccess() ) {
                        if ( tkr1ZDir == 0 )
                            tkr1ZDir = -1;
                    }

                    float bias = m_ubInterpolateTool->interpolate("Prof",log10(CAL_cfp_energy), tkr1ZDir);                  
                    CAL_cfp_energyUB = bias == 0 ? -1 : CAL_cfp_energy / bias;
                    // std::cout << "CalValsTool CAL_cfp_energy " << CAL_cfp_energy << " ( " << log10(CAL_cfp_energy) << " ) " << CAL_cfp_energyUB << " ( " << log10(CAL_cfp_energyUB) << " ) " << tkr1ZDir << std::endl;
                }
            }
            else if (corResult.getCorrectionName() == "NewCalFullProfileTool")
              {
                CAL_newcfp_energy    = corResult.getParams().getEnergy();
                CAL_newcfp_energy_err    = corResult["energy_err"];
                CAL_newcfp_totChiSq  = corResult["totchisq"];
                CAL_newcfp_seltotChiSq  = corResult["seltotchisq"];
                CAL_newcfp_ChiSqdist  = corResult["chisqdist"];
                CAL_newcfp_Chidist  = corResult["chidist"];
                CAL_newcfp_calEffRLn = corResult["cal_eff_RLn"];
                CAL_newcfp_tkrRLn = corResult["tkr_RLn"];
                CAL_newcfp_par0 = corResult["par0"];
                CAL_newcfp_par1 = corResult["par1"];
                CAL_newcfp_parcf = corResult["parcf"];
                CAL_newcfp_xcor = corResult["xcor"];
                CAL_newcfp_ycor = corResult["ycor"];
                CAL_newcfp_alpha = corResult["alpha"];
                CAL_newcfp_tmax = corResult["tmax"];
                CAL_newcfp_fiterrflg = corResult["fitflag"];
                CAL_newcfp_nlayersat = corResult["nlayersat"];
                CAL_newcfp_nxtalsat = corResult["nxtalsat"];
                CAL_newcfp_nxtalsel = corResult["nxtalsel"];
                CAL_newcfp_p0 = corResult["recp0"];
                CAL_newcfp_p1 = corResult["recp1"];
                CAL_newcfp_p2 = corResult["recp2"];
                CAL_newcfp_v0 = corResult["recv0"];
                CAL_newcfp_v1 = corResult["recv1"];
                CAL_newcfp_v2 = corResult["recv2"];
                //
                CAL_newcfp_calfit_energy    = corResult["calfit_fit_energy"];
                CAL_newcfp_calfit_energy_err    = corResult["calfit_energy_err"];
                CAL_newcfp_calfit_totChiSq  = corResult["calfit_totchisq"];
                CAL_newcfp_calfit_seltotChiSq  = corResult["calfit_seltotchisq"];
                CAL_newcfp_calfit_ChiSqdist  = corResult["calfit_chisqdist"];
                CAL_newcfp_calfit_Chidist  = corResult["calfit_chidist"];
                CAL_newcfp_calfit_calEffRLn = corResult["calfit_cal_eff_RLn"];
                CAL_newcfp_calfit_par0 = corResult["calfit_par0"];
                CAL_newcfp_calfit_par1 = corResult["calfit_par1"];
                CAL_newcfp_calfit_parcf = corResult["calfit_parcf"];
                CAL_newcfp_calfit_xcor = corResult["calfit_xcor"];
                CAL_newcfp_calfit_ycor = corResult["calfit_ycor"];
                CAL_newcfp_calfit_alpha = corResult["calfit_alpha"];
                CAL_newcfp_calfit_tmax = corResult["calfit_tmax"];
                CAL_newcfp_calfit_fiterrflg = corResult["calfit_fitflag"];
                CAL_newcfp_calfit_widening = corResult["calfit_widening"];
                CAL_newcfp_calfit_nxtalsel = corResult["calfit_nxtalsel"];
                //
                if ( CAL_newcfp_energy == 0 ) {
                  CAL_newcfp_energyUB = 0;
                }
                else {
                  CAL_newcfp_energyUB = CAL_newcfp_energy;
                }
                // UB for full profile with cal-fit, for Cal Only
                if ( CAL_newcfp_calfit_energy  == 0 ) {
                    CAL_newcfp_calfit_energyub = 0;
                }
                else {

                  // cs: we need Z Dir of first/best cluster.
                  // not sure why we fill Energy before clusters params, 
                  // is all this energy stuff referred to the first cluster? I hope so...
                  Event::CalCluster* cal1 = firstCluster;
                  float cal1ZDir = -1*cal1->getMomParams().getAxis().z();
                    
                  float bias = m_ubInterpolateTool->interpolate("NewProfCalOnly", log10(CAL_newcfp_calfit_energy),   cal1ZDir );
                  CAL_newcfp_calfit_energyub = bias == 0 ? -1 : CAL_newcfp_calfit_energy / bias;
                  //std::cout << "CalValsTool CAL_newcfp_calfit_energy " << CAL_newcfp_calfit_energy << " ( " << log10(CAL_newcfp_calfit_energy) << " ) " 
                  //          << tkr1ZDir <<" -> "<< CAL_newcfp_calfit_energyub << " ( bias: " << bias << " ) "  std::endl;
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
                // std::cout << "CalValsTool CAL_LkHd_energy " << CAL_LkHd_energy << ' ' << CAL_LkHd_energyUB << std::endl;
            }
        }
    }

    // //Make sure we have valid cluster data and some energy
    // if (!pCals) return sc;
    // if (pCals->empty()) return sc;

    // // Extract the uber information which is located at the end of the list
    // Event::CalCluster* calCluster = pCals->back();

    CAL_NumXtals_uber = uberCluster->getXtalsParams().getNumXtals();
    CAL_NumTruncXtals_uber = uberCluster->getXtalsParams().getNumTruncXtals();

    CAL_energy_uber = uberCluster->getMomParams().getEnergy();

    CAL_xEcntr_uber   = uberCluster->getPosition().x();
    CAL_yEcntr_uber   = uberCluster->getPosition().y();
    CAL_zEcntr_uber   = uberCluster->getPosition().z();
    CAL_xdir_uber     = uberCluster->getDirection().x();
    CAL_ydir_uber     = uberCluster->getDirection().y();
    CAL_zdir_uber     = uberCluster->getDirection().z();

    // Get pos and dir determined using only the transverse position information
    CAL_xEcntr2_uber  = uberCluster->getFitParams().getCentroid().x();
    CAL_yEcntr2_uber  = uberCluster->getFitParams().getCentroid().y();
    CAL_zEcntr2_uber  = uberCluster->getFitParams().getCentroid().z();
    CAL_nsaturated    = uberCluster->getNumSaturatedXtals();
    CAL_xdir2_uber    = uberCluster->getFitParams().getAxis().x();
    CAL_ydir2_uber    = uberCluster->getFitParams().getAxis().y();
    CAL_zdir2_uber    = uberCluster->getFitParams().getAxis().z();
    
    // Add CalTransRms for Uber cluster for development
    CAL_MomTransRms_uber = uberCluster->getMomParams().getTransRms();
    CAL_MomLongRms_uber  = uberCluster->getMomParams().getLongRms();
    CAL_ClassGamProb_uber = uberCluster->getClassParams().getProb("gam");
    CAL_ClassHadProb_uber = uberCluster->getClassParams().getProb("had");
    CAL_ClassGhostProb_uber = uberCluster->getClassParams().getProb("ghost");
    CAL_ClassMipProb_uber = uberCluster->getClassParams().getProb("mip");

    // Fill Uber2 (=Uber without second cluster): by default it is the same as Uber
    CAL_Uber2_NumXtals = uber2Cluster->getXtalsParams().getNumXtals();
    CAL_Uber2_RawEnergySum = uber2Cluster->getXtalsParams().getXtalRawEneSum();
    CAL_Uber2_XtalEneMax = uber2Cluster->getXtalsParams().getXtalEneMax();
    CAL_Uber2_MomXCntr = uber2Cluster->getMomParams().getCentroid().x();
    CAL_Uber2_MomYCntr = uber2Cluster->getMomParams().getCentroid().y();
    CAL_Uber2_MomZCntr = uber2Cluster->getMomParams().getCentroid().z();
    CAL_Uber2_MomXDir = uber2Cluster->getMomParams().getAxis().x();
    CAL_Uber2_MomYDir = uber2Cluster->getMomParams().getAxis().y();
    CAL_Uber2_MomZDir = uber2Cluster->getMomParams().getAxis().z();
    CAL_Uber2_MomNumCoreXtals = uber2Cluster->getMomParams().getNumCoreXtals();
    CAL_Uber2_MomTransRms = uber2Cluster->getMomParams().getTransRms();
    CAL_Uber2_MomLongRms = uber2Cluster->getMomParams().getLongRms();

    CAL_rest_energy   = 0.;
    CAL_rest_numXtals = 0.;

    Event::CalClusterVec::iterator calClusIter = rawClusterVec.begin();
    calClusIter++; // skip 1st cluster
    while(calClusIter != rawClusterVec.end())
      {
        Event::CalCluster* cluster = *calClusIter;
        CAL_rest_energy += cluster->getXtalsParams().getXtalCorrEneSum();
        CAL_rest_numXtals += cluster->getXtalsParams().getNumXtals();
        calClusIter++;
      }

    CAL_EnergyRaw  = firstCluster->getXtalsParams().getXtalCorrEneSum();
    if(CAL_EnergyRaw<1.0) return sc;

    for(int i = 0; i<m_nLayers; i++) CAL_eLayer[i] = (*firstCluster)[i].getEnergy();

    CAL_Trans_Rms = firstCluster->getMomParams().getTransRms();

    float logRLn = 0; 
    if ((CAL_LAT_RLn - CAL_Cntr_RLn) > 0.0) {
        logRLn = log(CAL_LAT_RLn - CAL_Cntr_RLn);
    }
    if (logRLn > 0.0) {
      CAL_Long_Rms  = firstCluster->getMomParams().getLongRms() / logRLn;
    }
    CAL_LRms_Asym = firstCluster->getMomParams().getLongRmsAsym();

    // Variables referring to the first cluster---added after the restructuring
    // of the CAL moments analysis (Luca Baldini, Dec. 26, 2010).
    // Basic variables.
    CAL_Clu1_NumXtals = firstCluster->getXtalsParams().getNumXtals();
    CAL_Clu1_NumTruncXtals = firstCluster->getXtalsParams().getNumTruncXtals();
    CAL_Clu1_NumSaturatedXtals = firstCluster->getXtalsParams().getNumSaturatedXtals();
    CAL_Clu1_RawEnergySum = firstCluster->getXtalsParams().getXtalRawEneSum();
    CAL_Clu1_CorrEnergySum = firstCluster->getXtalsParams().getXtalCorrEneSum();
    CAL_Clu1_XtalEneMax = firstCluster->getXtalsParams().getXtalEneMax();
    CAL_Clu1_XtalEneRms = firstCluster->getXtalsParams().getXtalEneRms();
    CAL_Clu1_XtalEneSkewness = firstCluster->getXtalsParams().getXtalEneSkewness();
    // Variables from the moments analysis.
    CAL_Clu1_MomXCntr = firstCluster->getMomParams().getCentroid().x();
    CAL_Clu1_MomYCntr = firstCluster->getMomParams().getCentroid().y();
    CAL_Clu1_MomZCntr = firstCluster->getMomParams().getCentroid().z();
    CAL_Clu1_MomXDir = firstCluster->getMomParams().getAxis().x();
    CAL_Clu1_MomYDir = firstCluster->getMomParams().getAxis().y();
    CAL_Clu1_MomZDir = firstCluster->getMomParams().getAxis().z();
    CAL_Clu1_MomNumIterations = firstCluster->getMomParams().getNumIterations();
    CAL_Clu1_MomRunFlag = 0.;
    if (firstCluster->checkStatusBit(Event::CalCluster::MOMENTS)) { CAL_Clu1_MomRunFlag = 1.; }
    CAL_Clu1_MomNumCoreXtals = firstCluster->getMomParams().getNumCoreXtals();
    CAL_Clu1_MomTransRms = firstCluster->getMomParams().getTransRms();
    CAL_Clu1_MomLongRms = firstCluster->getMomParams().getLongRms();
    CAL_Clu1_MomLongRmsAysm = firstCluster->getMomParams().getLongRmsAsym();
    CAL_Clu1_MomLongSkewness = firstCluster->getMomParams().getLongSkewness();
    CAL_Clu1_MomCoreEneFrac = firstCluster->getMomParams().getCoreEnergyFrac();
    CAL_Clu1_MomFullLenght = firstCluster->getMomParams().getFullLength();
    CAL_Clu1_MomdEdxAve = firstCluster->getMomParams().getdEdxAverage();
    CAL_Clu1_Mom_Cov_Cxx = firstCluster->getMomParams().getxPosxPos();
    CAL_Clu1_Mom_Cov_Cyy = firstCluster->getMomParams().getyPosyPos();
    CAL_Clu1_Mom_Cov_Czz = firstCluster->getMomParams().getzPoszPos();
    CAL_Clu1_Mom_Cov_Dxx = 0.;
    CAL_Clu1_Mom_Cov_Dyy = 0.;
    CAL_Clu1_Mom_Cov_Dzz = 0.;
    CAL_Clu1_Mom_Cov_Dxy = 0.;
    CAL_Clu1_Mom_Cov_Dxz = 0.;
    CAL_Clu1_Mom_Cov_Dyz = 0.;
    CAL_Clu1_Mom_Cov_Sxx = 0.;
    CAL_Clu1_Mom_Cov_Syy = 0.;
    CAL_Clu1_Mom_Cov_Sxy = 0.;

    if (firstCluster->checkStatusBit(Event::CalCluster::MOMENTS)) {
      if (CAL_Clu1_MomZDir>0.001){ // exclude horizontal directions
        CLHEP::HepMatrix CalCov = firstCluster->getMomParams().getMomErrsTkrRep();  
        CAL_Clu1_Mom_Cov_Sxx = CalCov(2,2);
        CAL_Clu1_Mom_Cov_Syy = CalCov(4,4);
        CAL_Clu1_Mom_Cov_Sxy = CalCov(2,4);
      }
      CLHEP::HepMatrix AxisCov = firstCluster->getMomParams().getAxisErrs();
      CAL_Clu1_Mom_Cov_Dxx = AxisCov(1,1);
      CAL_Clu1_Mom_Cov_Dyy = AxisCov(2,2);
      CAL_Clu1_Mom_Cov_Dzz = AxisCov(3,3);
      CAL_Clu1_Mom_Cov_Dxy = AxisCov(1,2);
      CAL_Clu1_Mom_Cov_Dxz = AxisCov(1,3);
      CAL_Clu1_Mom_Cov_Dyz = AxisCov(2,3);
    }
    // Variables from Philippe's fit.
    CAL_Clu1_FitXCntr = firstCluster->getFitParams().getCentroid().x();
    CAL_Clu1_FitYCntr = firstCluster->getFitParams().getCentroid().y();
    CAL_Clu1_FitZCntr = firstCluster->getFitParams().getCentroid().z();
    CAL_Clu1_FitXDir = firstCluster->getFitParams().getAxis().x();
    CAL_Clu1_FitYDir = firstCluster->getFitParams().getAxis().y();
    CAL_Clu1_FitZDir = firstCluster->getFitParams().getAxis().z();
    CAL_Clu1_FitNumLayers = firstCluster->getFitParams().getFitLayers();
    CAL_Clu1_FitChiSquare = firstCluster->getFitParams().getChiSquare();
    // Variables from the Minimum Spanning Tree clustering.
    CAL_Clu1_MstMinEdgeLen = firstCluster->getMSTreeParams().getMinEdgeLength();
    CAL_Clu1_MstMaxEdgeLen = firstCluster->getMSTreeParams().getMaxEdgeLength();
    CAL_Clu1_MstAveEdgeLen = firstCluster->getMSTreeParams().getMeanEdgeLength();
    CAL_Clu1_MstRmsEdgeLen = firstCluster->getMSTreeParams().getRmsEdgeLength();
    CAL_Clu1_MstAveTruncEdgeLen = firstCluster->getMSTreeParams().getMeanEdgeLengthTrunc();
    CAL_Clu1_MstRmsTruncEdgeLen = firstCluster->getMSTreeParams().getRmsEdgeLengthTrunc();
    // Variables from the cluster classification.
    CAL_Clu1_ClassGamProb = firstCluster->getClassParams().getProb("gam");
    CAL_Clu1_ClassHadProb = firstCluster->getClassParams().getProb("had");
    CAL_Clu1_ClassGhostProb = firstCluster->getClassParams().getProb("ghost");
    CAL_Clu1_ClassMipProb = firstCluster->getClassParams().getProb("mip");


    // SSD Veto stuff here:
    //if (firstCluster->checkStatusBit(Event::CalCluster::MOMENTS)){ // not sure what is the best condition to check
    if(firstCluster->getMomParams().getNumIterations()>0){
      CAL_Clu1_SSDVeto = CalSSDEvaluation(firstCluster); // cs: first cluster
      CAL_Clu1_SSDVetoNoHitFlag = (int)floor(m_SSDVetoNoHitFlag + 0.5);;
      CAL_Clu1_VetoPlaneCrossed = (int)floor(m_VetoPlaneCrossed + 0.5); 
      CAL_Clu1_VetoHitFound     = (int)floor(m_VetoHitFound+ 0.5);
      CAL_Clu1_VetoTrials       = (int)floor(m_VetoTrials + 0.5);
      CAL_Clu1_VetoUnknown      = (int)floor(m_VetoUnknown + 0.5);
      CAL_Clu1_VetoDeadPlane    = (int)floor(m_VetoDeadPlane + 0.5);
      CAL_Clu1_VetoTruncated    = (int)floor(m_VetoTruncated + 0.5); 
      CAL_Clu1_VetoTower        = (int)floor(m_VetoTower + 0.5);   
      CAL_Clu1_VetoGapCorner    = (int)floor(m_VetoGapCorner + 0.5);
      CAL_Clu1_VetoGapEdge      = (int)floor(m_VetoGapEdge + 0.5);   
      CAL_Clu1_VetoBadCluster   = (int)floor(m_VetoBadCluster + 0.5);
    }
    else { 
      CAL_Clu1_SSDVeto       = _badFloat;   
      CAL_Clu1_SSDVetoNoHitFlag = _badFloat;
      CAL_Clu1_VetoPlaneCrossed = _badFloat;
      CAL_Clu1_VetoHitFound     = _badFloat;
      CAL_Clu1_VetoTrials       = _badFloat;
      CAL_Clu1_VetoUnknown      = _badFloat;
      CAL_Clu1_VetoDeadPlane    = _badFloat;
      CAL_Clu1_VetoTruncated    = _badFloat;
      CAL_Clu1_VetoTower        = _badFloat;
      CAL_Clu1_VetoGapCorner    = _badFloat;
      CAL_Clu1_VetoGapEdge      = _badFloat;
      CAL_Clu1_VetoBadCluster   = _badFloat;
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

    std::vector<Event::CalXtalRecData*>::iterator jlog;
    {
        // Find Xtal with max. energy
        CAL_Num_Xtals = (float)firstClusterXtalList.size();
        for( jlog=firstClusterXtalList.begin(); jlog != firstClusterXtalList.end(); ++jlog) {
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
        no_xtals=firstClusterXtalList.size();
    }
    int no_xtals_trunc=firstCluster->getNumTruncXtals();
    CAL_Xtal_Ratio= (no_xtals>0) ? float(no_xtals_trunc)/no_xtals : 0;
    CAL_Num_Xtals_Trunc = float(no_xtals_trunc); 

    // No use in continuing if too little energy in CAL
    //if(CAL_EnergyRaw < 5.) return sc;
    // See JIRA LPATE-48

    Point  cal_pos  = firstCluster->getPosition();
    Vector cal_dir  = firstCluster->getDirection();

    Point cor_cal_pos = cal_pos;
    // Apply the centroid correction using the cal axis. If there is a track, the correction is applied after in the code using the track axis.
    if(cal_dir.z()!=0) cor_cal_pos = firstCluster->getCorPosition(cal_dir);

    CAL_xEcntr      = cal_pos.x();
    CAL_yEcntr      = cal_pos.y();
    CAL_zEcntr      = cal_pos.z();
    CAL_xdir        = cal_dir.x();
    CAL_ydir        = cal_dir.y();
    CAL_zdir        = cal_dir.z();

    CAL_Clu1_MomXCntrCor = cor_cal_pos.x();
    CAL_Clu1_MomYCntrCor = cor_cal_pos.y();
    CAL_Clu1_MomZCntrCor = cor_cal_pos.z();

    // Get pos and dir determined using only the transverse position information
    CAL_xEcntr2        = firstCluster->getFitParams().getCentroid().x();
    CAL_yEcntr2        = firstCluster->getFitParams().getCentroid().y();
    CAL_zEcntr2        = firstCluster->getFitParams().getCentroid().z();
    CAL_posdir_chisq   = firstCluster->getFitParams().getChiSquare();
    CAL_posdir_nlayers = firstCluster->getFitParams().getFitLayers();
    CAL_nsaturated     = firstCluster->getNumSaturatedXtals();
    CAL_xdir2          = firstCluster->getFitParams().getAxis().x();
    CAL_ydir2          = firstCluster->getFitParams().getAxis().y();
    CAL_zdir2          = firstCluster->getFitParams().getAxis().z();

    // Get the lower and upper limits for the CAL in the installed towers
    double deltaX = 0.5*(m_xNum*m_towerPitch - m_calXWidth);
    double deltaY = 0.5*(m_yNum*m_towerPitch - m_calYWidth);
    double calXLo = m_tkrGeom->getLATLimit(0,LOW)  + deltaX;
    double calXHi = m_tkrGeom->getLATLimit(0,HIGH) - deltaX;
    double calYLo = m_tkrGeom->getLATLimit(1,LOW)  + deltaY;
    double calYHi = m_tkrGeom->getLATLimit(1,HIGH) - deltaY;

    // define the distance from the cal centroid to the nearest gap
    // the second argument of activeDist (view) does not seem to be used
    int tmp_view;
    CAL_TwrEdgeCntr = activeDist(cal_pos, tmp_view);
    CAL_TwrEdgeCntrCor = activeDist(cor_cal_pos, tmp_view);
    
    // collect the CAL edge energy
    // the edge is larger of the width and length of the layer
    // we can do better if this is at all useful.
    {
        // do the Cal[X/Y]PosRmsLL here
        double eRmsLL = 0;
        double xLL = 0;
        double xSqLL = 0;
        double yLL = 0;
        double ySqLL = 0;
        int    nRmsLL = 0;
        bool   doLL = (m_nLayers>1 && m_nCsI>1); // true for Glast, false for EGRET

        std::vector<Event::CalXtalRecData*>::iterator jlog;
        for( jlog=firstClusterXtalList.begin(); jlog != firstClusterXtalList.end(); ++jlog) {
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
    if(firstCluster->getRmsLong() < .1 ) { // Trap no-calc. condition
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

    bool do1 = false;
    bool do2 = false;
    if(num_tracks > 0) { 
        // Get the first track
        pTrack1 = pTracks->begin();
        track_1 = *pTrack1;
        if ((track_1->getStatusBits() & Event::TkrTrack::COSMICRAY)==0) do1 = true; //RJ , LSR
    }
    if(do1) {

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

        cor_cal_pos = firstCluster->getCorPosition(t1);
        CAL_Track_DOCA_Cor = (float)track1.docaOfPoint(cor_cal_pos);
        CAL_Clu1_MomXCntrCor = cor_cal_pos.x();
        CAL_Clu1_MomYCntrCor = cor_cal_pos.y();
        CAL_Clu1_MomZCntrCor = cor_cal_pos.z();
        CAL_TwrEdgeCntrCor = activeDist(cor_cal_pos, tmp_view);

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
            if ((track_1->getStatusBits() & Event::TkrTrack::COSMICRAY)==0) do2 = true; //RJ , LSR
        }
        if(do2) {

            // Get the start and direction 
            x2 = track_1->getInitialPosition();
            t2 = track_1->getInitialDirection();

            Point x2Top = x2 + (pos0.z()-x2.z())/t2.z() * t2; 
            CAL_Track_Sep = (x2Top - pos0).mag();
        }

        // If vertexed - use first vertex
        //std::cout << "pVerts " << pVerts << std::endl;
        if(pVerts) {
            if(pVerts->size()>0) {
                Event::TkrVertex* vertex = *(pVerts->begin()); 
                x0 = vertex->getPosition();
                t0 = vertex->getDirection();
            }
        }

        // try Bill's new vars... 
        {
            //make a map of xtal energy by doca
            // we need this so that we can drop the largest entries for the truncated calc.
            std::multimap<double, double> docaMap;
            std::multimap<double, double>::iterator mapIter;

            std::vector<Event::CalXtalRecData*>::iterator jlog;
            for( jlog=firstClusterXtalList.begin(); jlog != firstClusterXtalList.end(); ++jlog) {
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
    
    if(num_tracks>0 && fabs(cal_dir.x()) < 1. &&
       firstCluster->checkStatusBit(Event::CalCluster::MOMENTS) ) {
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
    TSaxisP = firstCluster->getPosition();
    Vector TSaxisVin = firstCluster->getDirection();
    TSaxisV = TSaxisVin.unit();
    //
    // Filling TSdist...
    //
    TSfillTSdist(firstClusterXtalList);

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
        TSfillTSdist(firstClusterXtalList);

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
        TSfillTSdist(firstClusterXtalList);

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

    // for uber cluster
    if(mynumtracks>0 && uberCluster)
      { 
        m_pCalClusterHitTool->fillRecDataVec(uberCluster);
        std::vector<Event::CalXtalRecData*> uberClusterXtalList = m_pCalClusterHitTool->getRecDataVec();

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
        TSfillTSdist(uberClusterXtalList);
        
        if(TSnlog>0)
          {
            //
            // fill CAL_TS_TKR_T_
            //
            TSfillTS(0);
            CAL_uber_TS_TKR_T_100 = (float)TSTS[TSnlog-1];
            //
            // fill CAL_TS_TKR_TL_
            //
            TSfillTS(1);
            CAL_uber_TS_TKR_TL_100 = (float)TSTS[TSnlog-1];
        }
    }

    // for uber2 cluster
    if(mynumtracks>0 && uber2Cluster)
      { 
        m_pCalClusterHitTool->fillRecDataVec(uber2Cluster);
        std::vector<Event::CalXtalRecData*> uber2ClusterXtalList = m_pCalClusterHitTool->getRecDataVec();

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
        TSfillTSdist(uber2ClusterXtalList);
        
        if(TSnlog>0)
          {
            //
            // fill CAL_TS_TKR_T_
            //
            TSfillTS(0);
            CAL_uber2_TS_TKR_T_100 = (float)TSTS[TSnlog-1];
            //
            // fill CAL_TS_TKR_TL_
            //
            TSfillTS(1);
            CAL_uber2_TS_TKR_TL_100 = (float)TSTS[TSnlog-1];
        }
    }

    // Fill information of the clusters (on top of Clu1)

  CAL_Clu2_NumXtals = 0;
  CAL_Clu2_RawEnergySum = 0;
  CAL_Clu2_XtalEneMax = 0;
  CAL_Clu2_MomXCntr = 0;
  CAL_Clu2_MomYCntr = 0;
  CAL_Clu2_MomZCntr = 0;
  CAL_Clu2_MomXDir = 0;
  CAL_Clu2_MomYDir = 0;
  CAL_Clu2_MomZDir = 0;
  CAL_Clu2_MomNumCoreXtals = 0;
  CAL_Clu2_MomTransRms = 0;
  CAL_Clu2_MomLongRms = 0;
  CAL_Clu2_ClassGamProb = 0;
  CAL_Clu2_ClassHadProb = 0;
  CAL_Clu2_ClassGhostProb = 0;
  CAL_Clu2_ClassMipProb = 0;
  CAL_Clu2_Dist = 0;
  CAL_Clu2_Doca = 0;

  CAL_Clu3_NumXtals = 0;
  CAL_Clu3_RawEnergySum = 0;
  CAL_Clu3_XtalEneMax = 0;
  CAL_Clu3_MomXCntr = 0;
  CAL_Clu3_MomYCntr = 0;
  CAL_Clu3_MomZCntr = 0;
  CAL_Clu3_MomXDir = 0;
  CAL_Clu3_MomYDir = 0;
  CAL_Clu3_MomZDir = 0;
  CAL_Clu3_MomNumCoreXtals = 0;
  CAL_Clu3_MomTransRms = 0;
  CAL_Clu3_MomLongRms = 0;
  CAL_Clu3_ClassGamProb = 0;
  CAL_Clu3_ClassHadProb = 0;
  CAL_Clu3_ClassGhostProb = 0;
  CAL_Clu3_ClassMipProb = 0;
  CAL_Clu3_Dist = 0;
  CAL_Clu3_Doca = 0;

  CAL_Clu4_NumXtals = 0;
  CAL_Clu4_RawEnergySum = 0;
  CAL_Clu4_XtalEneMax = 0;
  CAL_Clu4_MomXCntr = 0;
  CAL_Clu4_MomYCntr = 0;
  CAL_Clu4_MomZCntr = 0;
  CAL_Clu4_MomXDir = 0;
  CAL_Clu4_MomYDir = 0;
  CAL_Clu4_MomZDir = 0;
  CAL_Clu4_MomNumCoreXtals = 0;
  CAL_Clu4_MomTransRms = 0;
  CAL_Clu4_MomLongRms = 0;
  CAL_Clu4_ClassGamProb = 0;
  CAL_Clu4_ClassHadProb = 0;
  CAL_Clu4_ClassGhostProb = 0;
  CAL_Clu4_ClassMipProb = 0;
  CAL_Clu4_Dist = 0;
  CAL_Clu4_Doca = 0;

  CAL_Clu5_NumXtals = 0;
  CAL_Clu5_RawEnergySum = 0;
  CAL_Clu5_XtalEneMax = 0;
  CAL_Clu5_MomXCntr = 0;
  CAL_Clu5_MomYCntr = 0;
  CAL_Clu5_MomZCntr = 0;
  CAL_Clu5_MomXDir = 0;
  CAL_Clu5_MomYDir = 0;
  CAL_Clu5_MomZDir = 0;
  CAL_Clu5_MomNumCoreXtals = 0;
  CAL_Clu5_MomTransRms = 0;
  CAL_Clu5_MomLongRms = 0;
  CAL_Clu5_ClassGamProb = 0;
  CAL_Clu5_ClassHadProb = 0;
  CAL_Clu5_ClassGhostProb = 0;
  CAL_Clu5_ClassMipProb = 0;
  CAL_Clu5_Dist = 0;
  CAL_Clu5_Doca = 0;

  int iclu = -1;
  int numClusters = 0;
  double lambda;

  calClusIter = rawClusterVec.begin();
  while(calClusIter != rawClusterVec.end())
    {
      Event::CalCluster* cluster = *calClusIter;
      calClusIter++;
      iclu++;
      //
      if(iclu==1)
        {
          CAL_Clu2_NumXtals = cluster->getXtalsParams().getNumXtals();
          CAL_Clu2_RawEnergySum = cluster->getXtalsParams().getXtalRawEneSum();
          CAL_Clu2_XtalEneMax = cluster->getXtalsParams().getXtalEneMax();
          CAL_Clu2_MomXCntr = cluster->getMomParams().getCentroid().x();
          CAL_Clu2_MomYCntr = cluster->getMomParams().getCentroid().y();
          CAL_Clu2_MomZCntr = cluster->getMomParams().getCentroid().z();
          CAL_Clu2_MomXDir = cluster->getMomParams().getAxis().x();
          CAL_Clu2_MomYDir = cluster->getMomParams().getAxis().y();
          CAL_Clu2_MomZDir = cluster->getMomParams().getAxis().z();
          CAL_Clu2_MomNumCoreXtals = cluster->getMomParams().getNumCoreXtals();
          CAL_Clu2_MomTransRms = cluster->getMomParams().getTransRms();
          CAL_Clu2_MomLongRms = cluster->getMomParams().getLongRms();
          CAL_Clu2_ClassGamProb = cluster->getClassParams().getProb("gam");
          CAL_Clu2_ClassHadProb = cluster->getClassParams().getProb("had");
          CAL_Clu2_ClassGhostProb = cluster->getClassParams().getProb("ghost");
          CAL_Clu2_ClassMipProb = cluster->getClassParams().getProb("mip");
          CAL_Clu2_Dist = sqrt((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*(cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)+
                               (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*(cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)+
                               (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*(cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr));
          lambda = ((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*CAL_Clu1_MomXDir+
                    (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*CAL_Clu1_MomYDir+
                    (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*CAL_Clu1_MomZDir);
          CAL_Clu2_Doca = sqrt((cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))*(cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))+
                               (cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))*(cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))+
                               (cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir))*(cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir)));
        }
      else if(iclu==2)
        {
          CAL_Clu3_NumXtals = cluster->getXtalsParams().getNumXtals();
          CAL_Clu3_RawEnergySum = cluster->getXtalsParams().getXtalRawEneSum();
          CAL_Clu3_XtalEneMax = cluster->getXtalsParams().getXtalEneMax();
          CAL_Clu3_MomXCntr = cluster->getMomParams().getCentroid().x();
          CAL_Clu3_MomYCntr = cluster->getMomParams().getCentroid().y();
          CAL_Clu3_MomZCntr = cluster->getMomParams().getCentroid().z();
          CAL_Clu3_MomXDir = cluster->getMomParams().getAxis().x();
          CAL_Clu3_MomYDir = cluster->getMomParams().getAxis().y();
          CAL_Clu3_MomZDir = cluster->getMomParams().getAxis().z();
          CAL_Clu3_MomNumCoreXtals = cluster->getMomParams().getNumCoreXtals();
          CAL_Clu3_MomTransRms = cluster->getMomParams().getTransRms();
          CAL_Clu3_MomLongRms = cluster->getMomParams().getLongRms();
          CAL_Clu3_ClassGamProb = cluster->getClassParams().getProb("gam");
          CAL_Clu3_ClassHadProb = cluster->getClassParams().getProb("had");
          CAL_Clu3_ClassGhostProb = cluster->getClassParams().getProb("ghost");
          CAL_Clu3_ClassMipProb = cluster->getClassParams().getProb("mip");
          CAL_Clu3_Dist = sqrt((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*(cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)+
                               (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*(cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)+
                               (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*(cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr));
          lambda = ((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*CAL_Clu1_MomXDir+
                    (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*CAL_Clu1_MomYDir+
                    (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*CAL_Clu1_MomZDir);
          CAL_Clu3_Doca = sqrt((cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))*(cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))+
                               (cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))*(cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))+
                               (cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir))*(cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir)));
        }
      else if(iclu==3)
        {
          CAL_Clu4_NumXtals = cluster->getXtalsParams().getNumXtals();
          CAL_Clu4_RawEnergySum = cluster->getXtalsParams().getXtalRawEneSum();
          CAL_Clu4_XtalEneMax = cluster->getXtalsParams().getXtalEneMax();
          CAL_Clu4_MomXCntr = cluster->getMomParams().getCentroid().x();
          CAL_Clu4_MomYCntr = cluster->getMomParams().getCentroid().y();
          CAL_Clu4_MomZCntr = cluster->getMomParams().getCentroid().z();
          CAL_Clu4_MomXDir = cluster->getMomParams().getAxis().x();
          CAL_Clu4_MomYDir = cluster->getMomParams().getAxis().y();
          CAL_Clu4_MomZDir = cluster->getMomParams().getAxis().z();
          CAL_Clu4_MomNumCoreXtals = cluster->getMomParams().getNumCoreXtals();
          CAL_Clu4_MomTransRms = cluster->getMomParams().getTransRms();
          CAL_Clu4_MomLongRms = cluster->getMomParams().getLongRms();
          CAL_Clu4_ClassGamProb = cluster->getClassParams().getProb("gam");
          CAL_Clu4_ClassHadProb = cluster->getClassParams().getProb("had");
          CAL_Clu4_ClassGhostProb = cluster->getClassParams().getProb("ghost");
          CAL_Clu4_ClassMipProb = cluster->getClassParams().getProb("mip");
          CAL_Clu4_Dist = sqrt((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*(cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)+
                               (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*(cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)+
                               (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*(cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr));
          lambda = ((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*CAL_Clu1_MomXDir+
                    (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*CAL_Clu1_MomYDir+
                    (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*CAL_Clu1_MomZDir);
          CAL_Clu4_Doca = sqrt((cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))*(cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))+
                               (cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))*(cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))+
                               (cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir))*(cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir)));
        }
      else if(iclu==4)
        {
          CAL_Clu5_NumXtals = cluster->getXtalsParams().getNumXtals();
          CAL_Clu5_RawEnergySum = cluster->getXtalsParams().getXtalRawEneSum();
          CAL_Clu5_XtalEneMax = cluster->getXtalsParams().getXtalEneMax();
          CAL_Clu5_MomXCntr = cluster->getMomParams().getCentroid().x();
          CAL_Clu5_MomYCntr = cluster->getMomParams().getCentroid().y();
          CAL_Clu5_MomZCntr = cluster->getMomParams().getCentroid().z();
          CAL_Clu5_MomXDir = cluster->getMomParams().getAxis().x();
          CAL_Clu5_MomYDir = cluster->getMomParams().getAxis().y();
          CAL_Clu5_MomZDir = cluster->getMomParams().getAxis().z();
          CAL_Clu5_MomNumCoreXtals = cluster->getMomParams().getNumCoreXtals();
          CAL_Clu5_MomTransRms = cluster->getMomParams().getTransRms();
          CAL_Clu5_MomLongRms = cluster->getMomParams().getLongRms();
          CAL_Clu5_ClassGamProb = cluster->getClassParams().getProb("gam");
          CAL_Clu5_ClassHadProb = cluster->getClassParams().getProb("had");
          CAL_Clu5_ClassGhostProb = cluster->getClassParams().getProb("ghost");
          CAL_Clu5_ClassMipProb = cluster->getClassParams().getProb("mip");
          CAL_Clu5_Dist = sqrt((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*(cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)+
                               (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*(cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)+
                               (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*(cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr));
          lambda = ((cluster->getMomParams().getCentroid().x()-CAL_Clu1_MomXCntr)*CAL_Clu1_MomXDir+
                    (cluster->getMomParams().getCentroid().y()-CAL_Clu1_MomYCntr)*CAL_Clu1_MomYDir+
                    (cluster->getMomParams().getCentroid().z()-CAL_Clu1_MomZCntr)*CAL_Clu1_MomZDir);
          CAL_Clu5_Doca = sqrt((cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))*(cluster->getMomParams().getCentroid().x()-(CAL_Clu1_MomXCntr+lambda*CAL_Clu1_MomXDir))+
                               (cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))*(cluster->getMomParams().getCentroid().y()-(CAL_Clu1_MomYCntr+lambda*CAL_Clu1_MomYDir))+
                               (cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir))*(cluster->getMomParams().getCentroid().z()-(CAL_Clu1_MomZCntr+lambda*CAL_Clu1_MomZDir)));
        }     
    }

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

int CalValsTool::TSfillTSdist(std::vector<Event::CalXtalRecData*> xtallist)
{
    int i;
    TSnlog = 0;
    for(i=0;i<1536;++i)
    {
        TSdistTL[i] = 0;
        TSdistT[i] = 0;
        TSenergy[i] = 0;
    }

    int itow,ilay,icol,itowx,itowy;
    double lambda;
    Point TSxtalP;
    Point TSxtalC;
    Vector TSxtalV;
    Vector TSTC;
    double lambdamax = 326./2;

    std::vector<Event::CalXtalRecData*>::iterator jlog;
    for( jlog=xtallist.begin(); jlog != xtallist.end(); ++jlog)
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
    SmartDataPtr<Event::CalClusterMap>
      pCalClusterMap(m_pEventSvc,EventModel::CalRecon::CalClusterMap); 
    if(pCalClusterMap)
      {
        Event::CalClusterVec rawClusterVec = (*pCalClusterMap).get(EventModel::CalRecon::CalRawClusterVec);
        Event::CalCluster* firstCluster = rawClusterVec.front();
        if(firstCluster)
          CAL_EnergyRaw  = firstCluster->getXtalsParams().getXtalCorrEneSum();
      }
}


float CalValsTool::CalSSDEvaluation(const Event::CalCluster* cluster)
{// Method to compute the number of SSD Vetoes for the given cal direction
  MsgStream log(msgSvc(), name());

  // init to bad
  m_VetoPlaneCrossed = _badFloat; 
  m_VetoTrials = _badFloat;
  m_SSDVeto = _badFloat;
  m_SSDVetoNoHitFlag = _badFloat;
  m_VetoUnknown = _badFloat;
  m_VetoDeadPlane = _badFloat;
  m_VetoTruncated = _badFloat; 
  m_VetoTower = -1.;   
  m_VetoGapCorner = _badFloat;
  m_VetoGapEdge = _badFloat;   
  m_VetoBadCluster = _badFloat;
  m_VetoHitFound = _badFloat;

  Point  xc = cluster->getMomParams().getCentroid(); // let start from cal centroid
  Vector t1 = -1*cluster->getMomParams().getAxis();  // notice -1: different convention for Tkr and Cal

  CLHEP::HepMatrix CalCov = cluster->getMomParams().getMomErrsTkrRep();
  // check on validity of the mom analysis calc on num core xtals:
  // if it is less than 3 it is not a valid fit and errors are meaningless
  if( cluster->getMomParams().getNumCoreXtals()<3 ) {return _badFloat;}

  // START from above the CAL (NO MCS in CsI)
  // Cov matrix should be propagated accordingly...
  double zRefAboveCal = 0.;
  Point  x1((zRefAboveCal-xc.z())*(t1.x()/t1.z()) + xc.x(), 
            (zRefAboveCal-xc.z())*(t1.y()/t1.z()) + xc.y(), 
            zRefAboveCal);


  // Check we start in side the LAT
  if(!m_tkrGeom->isInActiveLAT(x1) ) {return _badFloat;}

  // fill 4x4 representation of error matrix, these are the only non-zero elements.
  // Tkr rep. is (xPos, xSlp, yPos, ySlp)
  // Propagated to x1 - I hope this is right!
  
  double xPosxPos = (zRefAboveCal-xc.z())*CalCov(2,2) + CalCov(1,1);
  double yPosyPos = (zRefAboveCal-xc.z())*CalCov(4,4) + CalCov(3,3);
  double xSlpxSlp = CalCov(2,2);
  double ySlpySlp = CalCov(4,4);
  double xSlpySlp = CalCov(2,4);

  /*
  // Need to fake the Track Params object that will be used by the propagator
  // Eventually this will be the cal covariance matrix

  // Start with reasonable error, NO FULL Cov Matrix yet
  // Assuming error only on Theta (~psf) and x-y position, all uncorrelated
  double dtheta = m_vetoDTheta; // error on Theta ~ few deg - JO param
  double dx     = 0.;    // error on X ~ 1mm
  double dy     = 0.;    // error on Y ~ 1mm
  double dz     = 0.;    // error on Z ~ 1mm for projected x,y
  // (Theta,Phi) => t1=(xdir,ydir, zdir)
  // t1 + dx, dy => ((x0 , sx , y0 , sy))

  // fill 4x4 representation, these are the only non-zero elements.
  double xPosxPos = dx*dx;
  double yPosyPos = dy*dy;
  double xSlpxSlp = 0.;
  double ySlpySlp = 0.;
  double xSlpySlp = 0.;
  if (t1.z()!=1 && t1.z()!=0){ // protection for vertical&horizontal case (that should never happen)
  double CommonFactor = (1./(t1.z()*t1.z()*t1.z()*t1.z()) + 1./(t1.z()*t1.z()) + 1./(1. - t1.z()*t1.z()));
  xSlpxSlp = dtheta*dtheta*t1.x()*t1.x()*CommonFactor;
  ySlpySlp = dtheta*dtheta*t1.y()*t1.y()*CommonFactor;
  xSlpySlp = dtheta*dtheta*t1.x()*t1.y()*CommonFactor;
  } 
  */


  const Event::TkrTrackParams params(x1.x(), t1.x()/t1.z(), x1.y(), t1.y()/t1.z(),
                                     xPosxPos,     0.,     0.,       0,
                                                xSlpxSlp,  0.,     xSlpySlp,
                                                        yPosyPos,    0.,
                                                                   ySlpySlp);


  // Let's start with Sum of Raw Energy - check other possibilities
  float ConEne  = cluster->getXtalsParams().getXtalRawEneSum();
  
  int topPlane = m_tkrGeom->numPlanes()-1; 
  double topOfTkr = m_tkrGeom->getPlaneZ(topPlane) + 1.0;
  double arc_min = fabs((topOfTkr-x1.z())/t1.z());
  arc_min = std::min( arc_min, maxPath);

  bool upwards = true;

  try {
    m_G4PropTool->setStepStart(params, x1.z(), upwards);
    m_G4PropTool->step(arc_min);
  } catch( std::exception& /*e*/) {
    printHeader(log);
    setAnaTupBit();
    log << "See previous exception printout." << endreq;
    log << " Skipping the Cal1 Veto calculations" << endreq;
    return m_SSDVeto;

  } catch (...) {
    printHeader(log);
    setAnaTupBit();
    log << "Unknown exception, see previous exception message, if any" << endreq;
    log << "Skipping the CAL1 Veto calculations" << endreq;
    log << "Initial track parameters: pos: " << x1 << endreq 
        << "dir: " << t1 << " arclen: " << arc_min << endreq;
    return m_SSDVeto;
    }
    
  m_VetoPlaneCrossed = 0.0; 
  m_VetoTrials = 0.0;
  m_SSDVeto = 0.0;
  m_SSDVetoNoHitFlag = 0.0;
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
  
  // for the footprints
  double secthX = 1./sqrt(1.0 - t1.x()*t1.x());
  double secthY = 1./sqrt(1.0 - t1.y()*t1.y());
  
  double xVetoRgn = _vetoRegion*secthX;
  double yVetoRgn = _vetoRegion*secthY;

  double arcLen = m_G4PropTool->getStepArcLen(0);
  bool isFirstPlane = true;
  int  previousplane = -1;

  for(int istep = 1; istep < numSteps; ++istep) { 
    volId = m_G4PropTool->getStepVolumeId(istep);
    volId.prepend(prefix);
    Point x_step       = m_G4PropTool->getStepPosition(istep);
    arcLen += m_G4PropTool->getStepArcLen(istep);

    bool forward = false;
    const Event::TkrTrackParams newParams = 
      m_G4PropTool->getTrackParams(arcLen, ConEne, forward);

    // we're outside the LAT
    if(x_step.z() > topOfTkr || !m_tkrGeom->isInActiveLAT(x_step) ) break; 

    // check that it's really a TKR hit
    if(volId.size() != 9) continue; 
    if(!(volId[0]==0 && volId[3]==1)) continue; // !(LAT && TKR)
    if(volId[6]> 1) continue;  //It's a converter or some other tray element!
   
    //int tower = idents::TowerId(volId[2], volId[1]).id();
    int tray = volId[4];
    int view = volId[5];
    int face = volId[6];
    int layer = m_tkrGeom->trayToBiLayer(tray, face);
    int plane = m_tkrGeom->trayToPlane(tray, face);
    
    if(previousplane==plane) continue; // check we are not in the same plane as previous step
    previousplane = plane;
    m_VetoPlaneCrossed++;

   // now check if there's a hit near the extrapolated track!
   // if there is, reset the veto counter... we want leading non-hits

    int firstPlane = -1;
    if(isFirstPlane) {
      isFirstPlane = false;
      firstPlane = plane;
    }

    // from tkr code:
    
    m_VetoTrials = abs(plane-firstPlane) + 1;

    // Increase the search region if the propagated error is large enough
    double vetoRgn = (view==0 ? std::max(xVetoRgn, 2*sqrt(newParams(1,1))) : std::max(yVetoRgn, 2*sqrt(newParams(3,3))));

    int nVetoHits = pQueryClusters->numberOfHitsNear(view, layer, 
                                                     vetoRgn, x_step, t1);

    if (nVetoHits>0) {
      // found a hit, reset the SSDVeto
      m_SSDVeto = 0.0;
      m_SSDVetoNoHitFlag = 0.0;
      if(m_enableVetoDiagnostics) { m_VetoHitFound++; } 
    } else { 
      // no hit
      idents::TkrId tkrId(volId);
      Event::TkrTrackParams outParams;

      unsigned int status_bits = 0;
      int stage = -1;

      // look for reasons of mesing hists, gaps, bad hits, etc...
      stage = pFlagHits->flagHits(tkrId, newParams, x_step.z(),
                                  m_minVetoError, m_maxVetoError, m_vetoNSigma, 
                                  outParams, status_bits);

      m_SSDVetoNoHitFlag += 1.0;
      if(stage==-1) { 
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
