    
/** @file CalValsTool.cxx
@brief Calculates the Cal analysis variables
@author Bill Atwood, Leon Rochester

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/CalValsTool.cxx,v 1.65 2005/11/25 16:50:14 chamont Exp $
*/

// Include files


#include "ValBase.h"

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
#include "geometry/Ray.h"

#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagator.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "idents/TowerId.h" 
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/Rotation.h"

#include "geometry/Ray.h"

#include "TMath.h"

/*! @class CalValsTool
@brief calculates Cal values

  @authors Bill Atwood, Leon Rochester
  */
  
  class CalValsTool :   public ValBase
  {
  public:
      
      CalValsTool( const std::string& type, 
          const std::string& name, 
          const IInterface* parent);
      
      virtual ~CalValsTool() { }
      
      StatusCode initialize();
      
      StatusCode calculate();
      
  private:
      
      // some pointers to services  
      /// GlastDetSvc used for access to detector info
      IGlastDetSvc*    m_detSvc; 
      /// TkrGeometrySvc used for access to tracker geometry info
      ITkrGeometrySvc* m_tkrGeom;

      /// some Geometry
      double m_towerPitch;
      int    m_xNum;
      int    m_yNum;

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
      float CAL_No_Xtals_Trunc;
      float CAL_Long_Rms;
      float CAL_Trans_Rms;
      float CAL_LRms_Asym;

      float CAL_MIP_Diff; 
      float CAL_MIP_Ratio;

      // New variables for new energy correction tools
      float CAL_cfp_energy;      // Energy from Full Profile tool
      float CAL_cfp_totChiSq;    // Total ChiSquare of fit divided by 11
      float CAL_cfp_calEffRLn;   // Effective radiation lengths in the Cal

      float CAL_lll_energy;      // Energy from the Last Layer Likelihood tool
      float CAL_lll_energyErr;   // Chisquare from the Last Layer Likelihood
      float CAL_tkl_energy;      // Energy from the tracker Likelihood tool
      float CAL_tkl_energyErr;   // Chisquare from the tracker likelihood
      
      //Calimeter items with Recon - Tracks
      float CAL_Track_DOCA;
      float CAL_Track_Angle;
	  float CAL_Track_Sep;
      float CAL_x0;
      float CAL_y0;
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

          // find TkrGeometrySvc service
          if (service("TkrGeometrySvc", m_tkrGeom, true).isFailure()){
              log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
              return StatusCode::FAILURE;
          }

          m_towerPitch = m_tkrGeom->towerPitch();
          m_xNum = m_tkrGeom->numXTowers();
          m_yNum = m_tkrGeom->numYTowers();

          if (getCalInfo().isFailure()) {
               log << MSG::ERROR << "Couldn't initialize the CAL constants" << endreq;
             return StatusCode::FAILURE;
          }
      } else {
          return StatusCode::FAILURE;
      }

      // load up the map

      /** @page anatup_vars 
      @section calvalstool CalValsTool Variables
<table>
<tr><th> Variable </th><th> Description
<tr><td> CalEnergyRaw 
<td>        Sum of the raw energies in all the crystals.  
            Includes estimate of missed energy due to zero-supression.  
    This replaces the variable CalEnergySum.  NEW! 
<tr><td> CalEnergyCorr 
<td>        Cal Energy corrected layer-by-layer for edges and leakage.  
            This replaces the variable CalEneSumCorr. NEW! 
<tr><td> CalLeakCorr 
<td>        Leakage correction: this is the contained fraction of the total energy 
            after edge corrections.  
<tr><td> CalEdgeCorr 
<td>        Effective layer-by-layer edge correction mainly due to the gaps 
             between Cal modules; multiplicative 
<tr><td> CalTotalCorr 
<td>        Global total correction. Includes effect due to dead material; 
            multiplicative 
<tr><td> CalCsIRLn 
<td>        Total radiation lengths in crystals, integrated along the 
            event axis (line connecting the first hit in the tracker to the CAL energy centroid) 
<tr><td> CalTotRLn 
<td>        Total radiation lengths in the CAL, integrated along the event axis. 
<tr><td> CalCntRLn 
<td>        Radiation lengths integrated along the event axis, up to energy centroid 
<tr><td> CalLATRLn 
<td>        Total radiation lengths integrated along the event axis 
            (including the tracker). 
<tr><td> CalDeadTotRat 
<td>        Ratio of radiation lengths in dead material to CalTotRLn 
<tr><td> CalDeadCntRat 
<td>        Ratio of radiation lengths in dead material up to energy centroid, 
    to CalCntRat 
<tr><td> CalTPred 
<td>        Model-predicted energy centroid in radiation lengths 
<tr><td> CalDeltaT 
<td>        Difference between measured and predicted energy centroids 
<tr><td> CalTwrEdge 
<td>        Distance of the entry point of the best track from the tower boundary, 
            measured at the top of the CAL. 
<tr><td> CalLATEdge 
<td>        Closest distance of track 1, projected to the top of the CAL, 
            to the edge of the CAL layer, taking non-square shape into account. 
            This is essentially the old merit skirt variable. 
<tr><td> CalEdgeEnergy 
<td>        The sum of the raw energies in each crystal for which the energy centroid 
            is within _deltaEdge (currently 50 mm) of the outside edge of one of 
            the outside CAL modules.
            This is an attempt at a "anti-coincidence counter" for the CAL.
<tr><td> CalTwrEdgeCntr 
<td>        Distance of the energy centroid from the nearest tower boundary.  
<tr><td> CalGapFraction 
<td>        Approximate fraction of the shower volumn which falls in inter-tower gaps. 
<tr><td> CalTrackSep 
<td>        Distance between impact points of two best tracks at CAL front face; 
            zero if only one track 
<tr><td> CalTrackDoca 
<td>        Distance between the projected vertex (or track if only one track) 
    and the energy centroid, evaluated at the z of the centroid. 
<tr><td> CalTrackAngle 
<td>        Angle between "gamma" direction in the tracker and direction of the CAL "track" 
<tr><td> CalELayerN, N=0,7 
<td>        Energy deposited in layer N of the CAL 
<tr><td> CalLyr0Ratio 
<td>        Ratio of CalELayer0 to CalEnergyRaw 
<tr><td> CalLyr7Ratio 
<td>        Ratio of CalELayer7 to CalEnergyRaw 
<tr><td> CalBkHalfRatio 
<td>        Ratio of total energy in back half of CAL (layers 4-7) to 
            CalEnergyRaw 
<tr><td> CalXtalsTrunc 
<td>        Number of CAL Xtals with > %1 of CalEnergyRaw (see CalXtalRatio) 
<tr><td> CalXtalRatio 
<td>        Ratio of number of Xtals with energy > 1% of CalEnergyRaw to 
            total number of struck Xtals in the event. 
<tr><td> CalXtalMaxEne 
<td>        Maximum energy found in a single Xtal
<tr><td> CalLongRms 
<td>        rms of the average of the 1st and 3rd shower moments. 
            Indicates the length of the measured shower along the shower axis. 
<tr><td> CalLRmsAsym 
<td>        The asymetry of the 1st and 3rd shower moments.  
            This should be close to zero. Because of ordering of moments it is slightly ... (??)
<tr><td> CalTransRms 
<td>        rms of transverse position measurements.
<tr><td> CalMIPDiff 
<td>        Difference between measured energy and that expected 
            from a minimum-ionizing particle 
<tr><td> CalMIPRatio 
<td>        Ratio of measured energy to that expected from a 
            minimum-ionizing particle 
<tr><td> Cal[X/Y/Z]Ecentr 
<td>        Energy centroid in [x/y/z]
<tr><td> Cal[X/Y/Z]Dir 
<td>        [x/y/z] direction cosine of CAL "track" 
<tr><td> Cal[X/Y]0 
<td>        [x/y] position of CAL "track" 
</table>

      */

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
      addItem("CalLyr0Ratio",  &CAL_Lyr0_Ratio);
      addItem("CalLyr7Ratio",  &CAL_Lyr7_Ratio);
      addItem("CalBkHalfRatio",&CAL_BkHalf_Ratio);

      addItem("CalXtalsTrunc", &CAL_No_Xtals_Trunc);
      addItem("CalXtalRatio",  &CAL_Xtal_Ratio);
      addItem("CalXtalMaxEne", &CAL_Xtal_maxEne);

	  addItem("CalTransRms",   &CAL_Trans_Rms);
      addItem("CalLongRms",    &CAL_Long_Rms);
      addItem("CalLRmsAsym",   &CAL_LRms_Asym);

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

      addItem("CalCfpEnergy",  &CAL_cfp_energy);
      addItem("CalCfpChiSq",   &CAL_cfp_totChiSq);
      addItem("CalCfpEffRLn",  &CAL_cfp_calEffRLn);
      addItem("CalLllEnergy",  &CAL_lll_energy);
      addItem("CalLllEneErr",  &CAL_lll_energyErr);
      addItem("CalTklEnergy",  &CAL_tkl_energy);
      addItem("CalTklEneErr",  &CAL_tkl_energyErr);

      zeroVals();
      
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

    // Recover pointer to CalEventEnergy
    Event::CalEventEnergyCol * calEventEnergyCol = 
        SmartDataPtr<Event::CalEventEnergyCol>(m_pEventSvc,EventModel::CalRecon::CalEventEnergyCol) ;
    Event::CalEventEnergy * calEventEnergy = 0 ;
    if ((calEventEnergyCol!=0)&&(!calEventEnergyCol->empty()))
        calEventEnergy = calEventEnergyCol->front() ;

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
 
	            CAL_x0 = corResult["CalTopX0"];
	            CAL_y0 = corResult["CalTopY0"];
            }
            else if (corResult.getCorrectionName() == "CalFullProfileTool")
            {
                CAL_cfp_energy    = corResult.getParams().getEnergy();
                CAL_cfp_totChiSq  = corResult["totchisq"];
                CAL_cfp_calEffRLn = corResult["cal_eff_RLn"];
            }
            else if (corResult.getCorrectionName() == "CalLastLayerLikelihoodTool")
            {
                CAL_lll_energy    = corResult.getParams().getEnergy();
                CAL_lll_energyErr = corResult.getParams().getEnergyErr();
            }
            else if (corResult.getCorrectionName() == "CalTkrLikelihoodTool")
            {
                CAL_tkl_energy    = corResult.getParams().getEnergy();
                CAL_tkl_energyErr = corResult.getParams().getEnergyErr();
            }
        }
    }

    //Make sure we have valid cluster data
    if (!pCals) return sc;
    if (pCals->empty()) return sc;
    
    Event::CalCluster* calCluster = pCals->front();
    
    CAL_EnergyRaw  = calCluster->getCalParams().getEnergy();
    for(int i = 0; i<8; i++) CAL_eLayer[i] = (*calCluster)[i].getEnergy();

    CAL_Trans_Rms = sqrt(calCluster->getRmsTrans()/CAL_EnergyRaw);
    if ((CAL_LAT_RLn - CAL_Cntr_RLn) > 0.0) {
        CAL_Long_Rms  = sqrt(calCluster->getRmsLong()/CAL_EnergyRaw) / 
		            log(CAL_LAT_RLn - CAL_Cntr_RLn);
    } else {
        CAL_Long_Rms = 0;
    }
    CAL_LRms_Asym = calCluster->getRmsLongAsym();

    if(CAL_EnergyRaw>0.0) {
        CAL_Lyr0_Ratio  = CAL_eLayer[0]/CAL_EnergyRaw;
        CAL_Lyr7_Ratio  = CAL_eLayer[7]/CAL_EnergyRaw;
        CAL_BkHalf_Ratio = (CAL_eLayer[4]+CAL_eLayer[5]+
                            CAL_eLayer[6]+CAL_eLayer[7])/CAL_EnergyRaw;
    }
    
    int no_xtals=0;
    CAL_Xtal_maxEne = 0.; 
    Event::CalXtalRecCol::const_iterator jlog;
    if (pxtalrecs) {
        // Find Xtal with max. energy
        for( jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog) {
            const Event::CalXtalRecData& recLog = **jlog;    
            double eneLog = recLog.getEnergy();
            if(eneLog > CAL_Xtal_maxEne) CAL_Xtal_maxEne = eneLog; 
        }
 
	// Number of Xtals
        no_xtals=pxtalrecs->size();
    }
    int no_xtals_trunc=calCluster->getNumTruncXtals();
    CAL_Xtal_Ratio= (no_xtals>0) ? float(no_xtals_trunc)/no_xtals : 0;
    CAL_No_Xtals_Trunc = float(no_xtals_trunc); 

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
    
    // Get the lower and upper limits for the CAL in the installed towers
    double deltaX = 0.5*(m_xNum*m_towerPitch - m_calXWidth);
    double deltaY = 0.5*(m_yNum*m_towerPitch - m_calYWidth);
    double calXLo = m_tkrGeom->getLATLimit(0,LOW)  + deltaX;
    double calXHi = m_tkrGeom->getLATLimit(0,HIGH) - deltaX;
    double calYLo = m_tkrGeom->getLATLimit(1,LOW)  + deltaY;
    double calYHi = m_tkrGeom->getLATLimit(1,HIGH) - deltaY;

	// Event axis locations relative to towers and LAT
    Point pos0(CAL_x0, CAL_y0, m_calZTop);
	int view;
    CAL_TwrEdgeTop  = activeDist(pos0, view);
    CAL_TwrEdgeCntr = activeDist(cal_pos, view);
	// Find the distance closest to an edge
    double dX = std::max(calXLo-CAL_x0, CAL_x0-calXHi);
    double dY = std::max(calYLo-CAL_y0, CAL_y0-calYHi);
    CAL_LATEdge = -std::max(dX, dY);
    // collect the CAL edge energy
    // the edge is larger of the width and length of the layer
    // we can do better if this is at all useful.
    if (pxtalrecs) {
        for( jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog) {
            const Event::CalXtalRecData& recLog = **jlog;    
            Point pos = recLog.getPosition();
            double xPos = pos.x(); double yPos = pos.y();
            double minX, minY;
            minX = std::min(abs(xPos-calXLo),abs(xPos-calXHi));
            minY = std::min(abs(yPos-calYLo),abs(yPos-calYHi));
            // for later... no time to check this out now!
            /*
            idents::CalXtalId id = recLog.getPackedId();
            if (id.isX()) {
                double minX = 
                    std::min(abs(xPos-(calXLo+m_deltaSide)),abs(xPos-(calXHi-m_deltaSide)));
                double minY = std::min(abs(yPos-calYLo),abs(yPos-calYHi));
            } else {
                double minX = std::min(abs(xPos-calXLo),abs(xPos-calXHi));
                double minY = 
                    std::min(abs(yPos-(calYLo+m_deltaSide)),abs(yPos-(calYHi-m_deltaSide)));
            }
            */
            if ((minX<_deltaEdge) ||  (minY<_deltaEdge)) {
                double eneLog = recLog.getEnergy();
                CAL_EdgeEnergy += eneLog;
            }
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
	if(num_tracks > 0) { 
        // Get the first track
        Event::TkrTrackColConPtr pTrack1 = pTracks->begin();
        Event::TkrTrack* track_1 = *pTrack1;
    
        // Get the start and direction 
        x1 = track_1->getInitialPosition();
        t1 = track_1->getInitialDirection();
		x0 = x1; 
		t0 = t1;

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
        if(pVerts) {
            Event::TkrVertex* vertex = *(pVerts->begin()); 
            x0 = vertex->getPosition();
            t0 = vertex->getDirection();
		}
    }

	 // Find the distance for energy centroid to track axis
    Vector x_diff = x0 - cal_pos;
    double x_diff_sq = x_diff*x_diff;
    double x_diff_t0 = x_diff*t0;
    CAL_Track_DOCA = sqrt(std::max(0.0, x_diff_sq - x_diff_t0*x_diff_t0));

	// Find angle between Track and Cal. Moment axis
	// Note: the direction in Cal is opposite to tracking!
    if(fabs(cal_dir.x()) < 1.) {
        double cosCalt0 = std::min(1., -t0*cal_dir); 
        cosCalt0 = std::max(cosCalt0, -1.);  // just in case...
        CAL_Track_Angle = acos(cosCalt0);
    }
    else CAL_Track_Angle = -.1f; 

	// Energy compared to muon 
    CAL_MIP_Diff    = CAL_EnergyRaw- 12.07*CAL_CsI_RLn;
    const double minRadLen = 0.1;
    CAL_MIP_Ratio   = CAL_EnergyRaw/(12.07*std::max(CAL_CsI_RLn*1., minRadLen));

	// Ratios of rad. len. in material other then CsI
    CAL_DeadTot_Rat = m_radLen_Stuff/std::max(minRadLen, CAL_Tot_RLn*1.);
    CAL_DeadCnt_Rat = m_radLen_CntrStuff/std::max(minRadLen, CAL_Cntr_RLn*1.); 

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
