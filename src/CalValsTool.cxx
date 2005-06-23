    
/** @file CalValsTool.cxx
@brief Calculates the Cal analysis variables
@author Bill Atwood, Leon Rochester

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/CalValsTool.cxx,v 1.54 2005/06/17 05:11:16 atwood Exp $
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
      double m_calZTop;
      double m_calZBot;

	  // Internal Variables
	  double m_radLen_CsI;
      double m_radLen_Stuff;
      double m_radLen_CntrStuff;

      IPropagator * m_G4PropTool; 
      IValsTool* m_pTkrTool;

      double CAL_EnergyRaw;
	  double CAL_EnergyCorr;
      double CAL_Leak_Corr;
      double CAL_Edge_Corr;
	  double CAL_Total_Corr; 
   
      double CAL_CsI_RLn;
      double CAL_Tot_RLn;
      double CAL_Cntr_RLn; 
	  double CAL_LAT_RLn; 
      double CAL_DeadTot_Rat;
      double CAL_DeadCnt_Rat; 
 
      double CAL_t_Pred; 
      double CAL_deltaT;
      double CAL_xEcntr;
      double CAL_yEcntr;
      double CAL_zEcntr;
      double CAL_xdir;
      double CAL_ydir;
      double CAL_zdir;

	  double CAL_Gap_Fraction;  
      double CAL_TwrEdgeCntr;
	  double CAL_TwrEdgeTop;
      double CAL_LATEdge; 

      double CAL_Lyr0_Ratio;
      double CAL_Lyr7_Ratio;
      double CAL_BkHalf_Ratio;

      double CAL_Xtal_Ratio;
      double CAL_Xtal_maxEne; 
      double CAL_eLayer[8];
      double CAL_No_Xtals_Trunc;
      double CAL_Long_Rms;
      double CAL_Trans_Rms;
      double CAL_LRms_Asym;

      double CAL_MIP_Diff; 
      double CAL_MIP_Ratio;
      
      //Calimeter items with Recon - Tracks
      double CAL_Track_DOCA;
      double CAL_Track_Angle;
      double CAL_x0;
      double CAL_y0;
  };
  
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

          // pick up the propagator
 
          IToolSvc* toolSvc = 0;
          if(service("ToolSvc", toolSvc, true).isFailure()) {
              log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
              return StatusCode::FAILURE;
          }
          if(!toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)) {
              log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
              return StatusCode::FAILURE;
          }

          sc = toolSvc->retrieveTool("TkrValsTool", m_pTkrTool);
          if( sc.isFailure() ) {
              log << MSG::ERROR << "Unable to find tool: " "TkrValsTool" << endreq;
              return sc;
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
      addItem("CalTrackDoca",  &CAL_Track_DOCA);
      addItem("CalTrackAngle", &CAL_Track_Angle);

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

	// Recover pointer to CalEventEnergy info  
    Event::CalEventEnergy* calEventEnergy = 
                 SmartDataPtr<Event::CalEventEnergy>(m_pEventSvc, EventModel::CalRecon::CalEventEnergy);

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

            if (corResult.getCorrectionName() == "CalValsCorrTool"){
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

                break;
            }
        }
    }

    //Make sure we have valid cluster data
    if (!pCals) return sc;
    
    Event::CalCluster* calCluster = pCals->front();
    
    CAL_EnergyRaw  = calCluster->getCalParams().getEnergy();
    for(int i = 0; i<8; i++) CAL_eLayer[i] = (*calCluster)[i].getEnergy();

	CAL_Trans_Rms = sqrt(calCluster->getRmsTrans()/CAL_EnergyRaw);
	CAL_Long_Rms  = sqrt(calCluster->getRmsLong()/CAL_EnergyRaw) / 
		            log(CAL_LAT_RLn - CAL_Cntr_RLn);
	CAL_LRms_Asym = calCluster->getRmsLongAsym();

    if(CAL_EnergyRaw>0.0) {
        CAL_Lyr0_Ratio  = CAL_eLayer[0]/CAL_EnergyRaw;
        CAL_Lyr7_Ratio  = CAL_eLayer[7]/CAL_EnergyRaw;
        CAL_BkHalf_Ratio = (CAL_eLayer[4]+CAL_eLayer[5]+
                            CAL_eLayer[6]+CAL_eLayer[7])/CAL_EnergyRaw;
    }
    
    // Find Xtal with max. energy
    CAL_Xtal_maxEne = 0.; 
    for( Event::CalXtalRecCol::const_iterator jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog)
	{
        const Event::CalXtalRecData& recLog = **jlog;    
        double eneLog = recLog.getEnergy();
        if(eneLog > CAL_Xtal_maxEne) CAL_Xtal_maxEne = eneLog; 
    }
 
	// Number of Xtals
	int no_xtals=pxtalrecs->size();
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

	// Compare Track (or vertex) to Cal soltion
	Point x0  =  cal_pos;
	Vector t0 = -cal_dir;
	if(calCluster->getRmsLong() < .1 ) { // Trap no-calc. condition
		x0 = Point(0., 0., 0.);
		t0 = Vector(0., 0., -1.);
	}
    
    // If there are tracks, use the best one
	int num_tracks = 0;
	if(pTracks) num_tracks = pTracks->size();
	if(num_tracks > 0) { 
    
        // Get the first track
        Event::TkrTrackColConPtr pTrack1 = pTracks->begin();
        Event::TkrTrack* track_1 = *pTrack1;
    
        // Get the start and direction 
        x0 = track_1->getInitialPosition();
        t0 = track_1->getInitialDirection();
    
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
    CAL_Track_DOCA = sqrt(x_diff_sq - x_diff_t0*x_diff_t0);

	// Find angle between Track and Cal. Moment axis
	// Note: the direction in Cal is opposite to tracking!
    if(fabs(cal_dir.x()) < 1.) {
        double cosCalt0 = -t0*cal_dir; 
        CAL_Track_Angle = acos(cosCalt0);
    }
    else CAL_Track_Angle = -.1; 

	// Energy compared to muon 
    CAL_MIP_Diff    = CAL_EnergyRaw- 12.07*CAL_CsI_RLn;
    const double minRadLen = 0.1;
    CAL_MIP_Ratio   = CAL_EnergyRaw/(12.07*std::max(CAL_CsI_RLn, minRadLen));

	// Ratios of rad. len. in material other then CsI
    CAL_DeadTot_Rat = m_radLen_Stuff/std::max(minRadLen, CAL_Tot_RLn);
    CAL_DeadCnt_Rat = m_radLen_CntrStuff/std::max(minRadLen, CAL_Cntr_RLn); 

    return sc;
}

StatusCode CalValsTool::getCalInfo()
{
    m_calZTop = m_tkrGeom->calZTop();
    m_calZBot = m_tkrGeom->calZBot();
    m_calXWidth = m_tkrGeom->calXWidth();
    m_calYWidth = m_tkrGeom->calYWidth();

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