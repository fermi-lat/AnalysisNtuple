/** @file CalValsTool.cxx
@brief Calculates the Cal analysis variables
@author Bill Atwood, Leon Rochester

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/CalValsTool.cxx,v 1.28 2003/08/24 22:39:47 atwood Exp $
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
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

#include "GlastSvc/Reco/IKalmanParticle.h"
#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagator.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "idents/TowerId.h" 
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "TMath.h"

#ifndef M_PI
#define M_PI = 3.14159265358979323846
#endif

namespace {
    
    double signC(double x) { return (x<0.) ? -1.:1.;} 
 
	double activeDist(double x, double y, double xpitch, double ypitch, int nx, int ny, int &XY, int &outer) {
        double edge = 0.; 
        double x_twr = signC(x)*(fmod(fabs(x),xpitch) - xpitch/2.);
        double y_twr = signC(y)*(fmod(fabs(y),ypitch) - ypitch/2.);
        
        outer = 0; 
        
        if(fabs(x_twr)/xpitch > fabs(y_twr)/ypitch) {
            edge = xpitch/2. - fabs(x_twr);
            XY = 1; 
            if(fabs(x) > 0.5*(nx-1)*xpitch) outer = 1;
        }
        else {
            edge = ypitch/2. - fabs(y_twr);
            XY = 2;
            if(fabs(y) > 0.5*(ny-1)*ypitch) outer = 1;
        }
        return edge;
    }
    
    double circle_frac(double r) {
        double rl = (fabs(r) < 1.) ? fabs(r):1.; 
        double a_slice = 2.*(M_PI/4. - rl*sqrt(1.-rl*rl)/2. - asin(rl)/2.);
        double in_frac = 1.-a_slice/M_PI;
        if(r < 0.) in_frac = a_slice/M_PI;
        return in_frac;
    }
    
    double circle_frac_simp(double r, double angle_factor) {
        double slice_0 = circle_frac(r);
        double slice_p = circle_frac(r+angle_factor);
        double slice_m = circle_frac(r-angle_factor);
        return (slice_p + 4.*slice_0 + slice_m)/6.;
    }
    
    double contained_frac(double x, double y, int nx, int ny, double pitch, double gap,  
        double r, double costh, double phi) {
        // Get the projected angles for the gap
        double tanth = sqrt(1.-costh*costh)/costh;
        double gap_x = gap - 40.*sin(phi)*tanth; //was 20 in RR185
        if(gap_x < 5.) gap_x = 5.;
        double gap_y = gap - 40.*cos(phi)*tanth;
        if(gap_y < 5.) gap_y = 5.; 
        
        // X Edges
        double x_twr = signC(x)*(fmod(fabs(x),pitch) - pitch/2.);
        double edge = pitch/2. - fabs(x_twr);
        double r_frac_plus = (edge-gap_x/2.)/r; 
        double angle_factor = sin(phi)*(1./costh - 1.);
        double in_frac_x  =  circle_frac_simp(r_frac_plus, angle_factor);
        if(fabs(x) < 0.5*(nx-1)*pitch) { // X edge is not outside limit of LAT
            double r_frac_minus = (edge + gap_x/2.)/r;
            in_frac_x += circle_frac_simp(-r_frac_minus, angle_factor);
        }
        
        // Y Edges
        double y_twr = signC(y)*(fmod(fabs(y),pitch) - pitch/2.);
        edge = pitch/2. - fabs(y_twr);
        r_frac_plus = (edge-gap_y/2.)/r; 
        angle_factor = cos(phi)*(1./costh - 1.);
        double in_frac_y  =  circle_frac_simp(r_frac_plus, angle_factor);
        if(fabs(y) < 0.5*(ny-1)*pitch) { // X edge is not outside limit of LAT
            double r_frac_minus = (edge + gap_y/2.)/r;
            in_frac_y += circle_frac_simp(-r_frac_minus, angle_factor);
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
    
    double beta(double E) {
        // Following expression for b (beta) comes from CalRecon
        // double beta      = exp(0.031*log(CAL_EnergySum/1000.))/2.29;
        // The above expression gives too small values
        double b      = 1.2*exp(0.031*log(E/1000.))/2.3;   //July 16, 2003 * 1.2
		//double b      = .51*exp(0.02*log(E/1000.));  // 17-July
        return b;
    }
    
    double alpha(double E) {
        // Following expression for a (alpha) comes from CalRecon
        // double alpha      = 2.65*exp(0.15*log(E/1000.));
        // The above expression gives too large values
        //double a = 2.50*exp(0.19*log(E/1000.));    //RR182
		// Found that these power law models diverge as E increases.  
		//double a = .98*log10(E) - .23; 
		double a = .9*log10(E); 
        return a;
	}
    
    double erf_cal(double x) {
        double t = 1./(1.+.47047*x);
        double results = 1. -(.34802 - (.09587 -.74785*t)*t)*t*exp(-x*x);
        return results;
    }
    double cal_trans(double x) {
        if(x < 0) return (.5*(1. + erf_cal(-x)));
        else      return (.5*(1. - erf_cal( x)));
    }
    
}

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
      ITkrGeometrySvc* m_geoSvc;

      /// some Geometry
      double m_towerPitch;
      int    m_xNum;
      int    m_yNum;
      /// gets the CAL info from detModel
      StatusCode getCalInfo();

      /// CAL vars
      double m_calXWidth;
      double m_calYWidth;
      double m_calZTop;
      double m_calZBot;

      IPropagatorSvc* m_propSvc;
	  IPropagator * m_G4PropTool; 


      
      //Global Calorimeter Tuple Items
      double CAL_EnergySum;
      //double CAL_Leak_Corr; 
      double CAL_Leak_Corr2;
      //double CAL_Edge_Corr; 
      double CAL_EdgeSum_Corr;     
      //double CAL_Total_Corr; 
      double CAL_TotSum_Corr; 
	  //double CAL_Energy_LLCorr; 

	  double CAL_CsI_RLn;
      double CAL_Tot_RLn;
      double CAL_Cnt_RLn; 
      double CAL_DeadTot_Rat;
      double CAL_DeadCnt_Rat; 
      //double CAL_a_Parm;
      //double CAL_b_Parm; 
      double CAL_t_Pred; 
      double CAL_deltaT;
      
      double CAL_EneSum_Corr;
      double CAL_Energy_Corr;
      double CAL_xEcntr;
      double CAL_yEcntr;
      double CAL_zEcntr;
      double CAL_xdir;
      double CAL_ydir;
      double CAL_zdir;

      double CAL_TwrEdge;
	  double CAL_LATEdge; 

      double CAL_TE_Nrm;
      double CAL_Track_Sep;

	  double CAL_Lyr0_Ratio;
	  double CAL_Lyr7_Ratio;
	  double CAL_BkHalf_Ratio;

	  double CAL_Xtal_Ratio;
	  double CAL_eLayer[8];
	  double CAL_No_Xtals_Trunc;
	  double CAL_Long_Rms;
	  double CAL_Trans_Rms;
	  double CAL_LRms_Ratio;

	  double CAL_MIP_Diff; 

      
      //Calimeter items with Recon - Tracks
      double CAL_Track_DOCA;
      double CAL_Track_Angle;
	  double CAL_TwrGap; 
      double CAL_x0;
      double CAL_y0;
      double CAL_z0;
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

          if( m_detSvc->getNumericConstByName("towerPitch", &m_towerPitch).isFailure() ||
              m_detSvc->getNumericConstByName("xNum",  &m_xNum).isFailure() ||
              m_detSvc->getNumericConstByName("yNum",  &m_yNum).isFailure() )  {
              log << MSG::ERROR << "Couldn't get detModel consts" << endreq;
              return StatusCode::FAILURE;
          }

          // find TkrGeometrySvc service
          if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
              log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
              return StatusCode::FAILURE;
          }

          if (getCalInfo().isFailure()) {
               log << MSG::ERROR << "Couldn't initialize the CAL constants" << endreq;
             return StatusCode::FAILURE;
          }

          // pick up the chosen propagator
          if (service("GlastPropagatorSvc", m_propSvc, true).isFailure()) {
              log << MSG::ERROR << "Couldn't find the GlastPropagatorSvc!" << endreq;
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
          return StatusCode::FAILURE;
      }
      
      // load up the map
      
      addItem("CalEnergySum",  &CAL_EnergySum);
      addItem("CalEnergyCorr", &CAL_Energy_Corr);
      addItem("CalEneSumCorr", &CAL_EneSum_Corr);
	  //addItem("Cal_Energy_LLCorr", &CAL_Energy_LLCorr);
   
      //addItem("CalLeakCorr",   &CAL_Leak_Corr);
      addItem("CalLeakCorr2",  &CAL_Leak_Corr2);
      
      //addItem("CalEdgeCorr",   &CAL_Edge_Corr);
      addItem("CalEdgeSumCorr",&CAL_EdgeSum_Corr);
      //addItem("CalTotalCorr",  &CAL_Total_Corr);
      addItem("CalTotSumCorr", &CAL_TotSum_Corr);
 
	  addItem("CalCsIRLn",     &CAL_CsI_RLn);
      addItem("CalTotRLn",     &CAL_Tot_RLn);
      addItem("CalCntRLn",     &CAL_Cnt_RLn);
      addItem("CalDeadTotRat", &CAL_DeadTot_Rat);
      addItem("CalDeadCntRat", &CAL_DeadCnt_Rat);
      //addItem("CalAParm",      &CAL_a_Parm);
      //addItem("CalBParm",      &CAL_b_Parm);
      addItem("CalTPred",      &CAL_t_Pred);
      addItem("CalDeltaT",     &CAL_deltaT);
      
      addItem("CalTwrEdge",    &CAL_TwrEdge);
	  addItem("CalLATEdge",    &CAL_LATEdge);
      addItem("CalTENrm",      &CAL_TE_Nrm);
      addItem("CalTrackSep",   &CAL_Track_Sep);
      addItem("CalTrackDoca",  &CAL_Track_DOCA);
	  addItem("CalTwrGap",     &CAL_TwrGap);

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
	  addItem("CalBkHalfRatio",  &CAL_BkHalf_Ratio);

	  addItem("CalXtalsTrunc", &CAL_No_Xtals_Trunc);
	  addItem("CalXtalRatio",  &CAL_Xtal_Ratio);

	  addItem("CalLongRms",    &CAL_Long_Rms);
	  addItem("CalLRmsRatio",  &CAL_LRms_Ratio);
	  addItem("CalTransRms",   &CAL_Trans_Rms);

	  addItem("CalMIPDiff",   &CAL_MIP_Diff);

      addItem("CalXEcntr",     &CAL_xEcntr);
      addItem("CalYEcntr",     &CAL_yEcntr);
      addItem("CalZEcntr",     &CAL_zEcntr);
      addItem("CalXDir",       &CAL_xdir);
      addItem("CalYDir",       &CAL_ydir);
      addItem("CalZDir",       &CAL_zdir);
      
      addItem("CalX0",         &CAL_x0);
      addItem("CalY0",         &CAL_y0);
      addItem("CalZ0",         &CAL_z0);

      
      zeroVals();
      
      return sc;
}


StatusCode CalValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    //MsgStream logstream(msgSvc(), name());

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>  
        pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    SmartDataPtr<Event::TkrVertexCol>     
        pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    //SmartDataPtr<Event::TkrClusterCol> 
    //    pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    SmartDataPtr<Event::CalClusterCol>     
        pCals(m_pEventSvc,EventModel::CalRecon::CalClusterCol);
	SmartDataPtr<Event::CalXtalRecCol> 
		pxtalrecs(m_pEventSvc,EventModel::CalRecon::CalXtalRecCol);
        
    //Make sure we have valid cluster data
    if (!pCals) return sc;
    
    Event::CalCluster* calCluster = pCals->front();
    
    CAL_EnergySum   = calCluster->getEnergySum();
	for(int i = 0; i<8; i++) CAL_eLayer[i] = calCluster->getEneLayer(i);
	CAL_Lyr0_Ratio  = CAL_eLayer[0]/CAL_EnergySum;
	CAL_Lyr7_Ratio  = CAL_eLayer[7]/CAL_EnergySum;
	CAL_BkHalf_Ratio = (CAL_eLayer[4]+CAL_eLayer[5]+CAL_eLayer[6]+CAL_eLayer[7])/CAL_EnergySum;

	CAL_Long_Rms      = calCluster->getRmsLong();
	CAL_Trans_Rms     = calCluster->getRmsTrans();
	CAL_LRms_Ratio    = CAL_Long_Rms / CAL_EnergySum;
	
	//CAL_Energy_LLCorr = calCluster->getEnergyCorrected();

	//Code from meritAlg
	int no_xtals=0;
    int no_xtals_trunc=0;
    for( Event::CalXtalRecCol::const_iterator jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog){

        const Event::CalXtalRecData& recLog = **jlog;
        
        double eneLog = recLog.getEnergy();
        if(eneLog>0)no_xtals++;
        if(eneLog>0.01*CAL_EnergySum)no_xtals_trunc++;
    }
    CAL_Xtal_Ratio= (no_xtals>0) ? float(no_xtals_trunc)/no_xtals : 0;
	CAL_No_Xtals_Trunc = float(no_xtals_trunc); 

    //Overwritten below, if there are tracks
    CAL_Energy_Corr = calCluster->getEnergyCorrected(); 
    if(CAL_EnergySum < 2.) return sc;  
    
    Point  cal_pos  = calCluster->getPosition();
    Vector cal_dir  = calCluster->getDirection();
    
    CAL_xEcntr       = cal_pos.x();
    CAL_yEcntr       = cal_pos.y();
    CAL_zEcntr       = cal_pos.z();
    CAL_xdir        = cal_dir.x();
    CAL_ydir        = cal_dir.y();
    CAL_zdir        = cal_dir.z();
    
    double twr_pitch = m_towerPitch;
    int iView = 0;
    int outside = 0;
    CAL_TwrEdge = activeDist(CAL_xEcntr, CAL_yEcntr, twr_pitch, twr_pitch, 
		                     m_xNum, m_yNum, iView, outside);
    
    if(iView==1) CAL_TE_Nrm  = cal_dir.x();
    else         CAL_TE_Nrm  = cal_dir.y();
    
    // with no tracks we've done what we can do.
    if(!pTracks) return sc; 
    int num_tracks = pTracks->size(); 
    if(num_tracks <= 0 ) return sc;
    
    // Get the first track
    Event::TkrFitConPtr pTrack1 = pTracks->begin();
    Event::TkrKalFitTrack* track_1 = dynamic_cast<Event::TkrKalFitTrack*>(*pTrack1);
    
    // Get the start and direction 
    Point  x0 = track_1->getPosition();
    Vector t0 = track_1->getDirection();
    
    // If vertexed - use first vertex
    if(pVerts) {
        Event::TkrVertexColPtr gammaPtr =  pVerts->begin(); 
        Event::TkrVertex *gamma = *gammaPtr; 
        x0 = gamma->getPosition();
        t0 = gamma->getDirection();
    }
    
    double costh  = fabs(t0.z()); 
    double phi    = atan(-t0.y()/t0.x());
    double phi_90 = fabs(phi); 
    
    double phi_45 = fabs(t0.y()/t0.x());
    if(phi_45 > 1.) phi_45 = 1./phi_45;
    phi_45 = atan(phi_45); 
    
    Vector x_diff = x0 - cal_pos;
    double x_diff_sq = x_diff*x_diff;
    double x_diff_t0 = x_diff*t0;
    CAL_Track_DOCA = sqrt(x_diff_sq - x_diff_t0*x_diff_t0);
    // The direction in Cal is opposite to tracking!
    double cosCalt0 = -t0*cal_dir; 
    CAL_Track_Angle = acos(cosCalt0);
    
    // Section to apply edge and leakage correction to Cal data. 
    // First apply layer by layer edge correction  
    Ray axis(x0, -t0); //Ray axis(cal_pos, x_diff.unit());
    double arc_len = -(m_calZTop-x0.z())/t0.z(); 
    Point z_0 = axis.position(arc_len); // Track entry point to top of Cal Stack 
    axis = Ray(z_0, x_diff.unit()); //Ray(z_0, -t0); //

    std::vector<Vector> pos_Layer = calCluster->getPosLayer();
    std::vector<double> ene_Layer = calCluster->getEneLayer();
    
    double ene_sum_corr = 0.;
    
    // Factors controlling edge correction: core radius, fringe radius, core fraction
    //    and size of inter-tower gap (active-to-active area) 
    double rm_hard   = 18.; // Fixed at ~ 1 rad. len in CsI 
    
    double rm_soft   = 55.+ 200.*cal_trans((CAL_EnergySum-250)/200.); 
     
    double hard_frac =.50 +.35*cal_trans((670.- CAL_EnergySum)/300.);   
      
    // is this related to towerPitch - CalModule width??
    double gap       = 48;
    
    // Reset shower area circle when multiple tracks are present
    // This increasingly matters below 1000 MeV.
	Event::TkrFitTrackBase::TrackEnd end = Event::TkrFitTrackBase::End;
	Ray trj_1 = Ray(track_1->getPosition(end), track_1->getDirection(end));
	double delta_z = trj_1.position().z() - m_calZTop;
    arc_len = delta_z/fabs(trj_1.direction().z());
    Point cal_1 = trj_1.position(arc_len); 
	if(fabs(cal_1.x()) > fabs(cal_1.y())) CAL_LATEdge = m_calXWidth/2. - fabs(cal_1.x());
    else                                  CAL_LATEdge = m_calYWidth/2. - fabs(cal_1.y());

    if(num_tracks > 1) {
        pTrack1++;
        Event::TkrKalFitTrack* track_2  
            = dynamic_cast<Event::TkrKalFitTrack*>(*pTrack1);
 
 
        Ray trj_2 = Ray(track_2->getPosition(end), track_2->getDirection(end));
 
        delta_z = trj_2.position().z() - m_calZTop;
        arc_len = delta_z/fabs(trj_2.direction().z());
        Point cal_2 = trj_2.position(arc_len);
        double separation = (cal_1.x()-cal_2.x())*(cal_1.x()-cal_2.x()) +
            (cal_1.y()-cal_2.y())*(cal_1.y()-cal_2.y());
        separation = sqrt(separation); 
        CAL_Track_Sep = separation;
        
        // Broaden shower core 
        if(separation/2. > rm_hard) { 
            rm_hard = separation/2.;
        }
        // Use track separation to limit the core fraction 
        double saturation = separation/75.; 
        if((1.-saturation) > hard_frac) hard_frac = 1. - saturation;
    }
    // Limit the hard(core) fraction
    if(hard_frac < .50) hard_frac = .50; 
    if(hard_frac > .85) hard_frac = .85; 
    
    double max_corr = 1.4;// + 3.6*cal_trans((80-CAL_EnergySum)/30.);
    
    // Apply circle correction layer by layer in the calorimeter
    // This is a essentially just a geometric correction 
    double edge_corr = 0.; 
    double good_layers = 0.; 
    for(int i=0; i<8; i++){
        if(ene_Layer[i] < 2.) {
            ene_sum_corr += ene_Layer[i];
            continue; 
        }
        Vector pos = z_0;
        double arc_len = (pos_Layer[i] - pos).magnitude();
        Point xyz_layer = axis.position(-arc_len);
        
        double in_frac_soft = contained_frac(xyz_layer.x(), xyz_layer.y(), 
            m_xNum, m_yNum, twr_pitch, gap, rm_soft, costh, phi_90);
        
        // Cut off correction upon leaving (through a side)
        if(in_frac_soft < .5) in_frac_soft = .5;
        
        double in_frac_hard = contained_frac(xyz_layer.x(), xyz_layer.y(), 
            m_xNum, m_yNum, twr_pitch, gap, rm_hard, costh, phi_90);
        
        double corr_factor 
            = 1./((1.-hard_frac)*in_frac_soft + hard_frac*in_frac_hard);
        if(corr_factor > max_corr) corr_factor = max_corr;  
        double ene_corr = ene_Layer[i]*corr_factor;
        ene_sum_corr += ene_corr;
        edge_corr    += corr_factor;
        good_layers  += 1.; 
    }
    if(ene_sum_corr < 1.) return sc;
    
    CAL_EdgeSum_Corr = ene_sum_corr/CAL_EnergySum;

    if (good_layers>0) edge_corr /= good_layers;
   // CAL_Edge_Corr = edge_corr; 
    
    // Set some Cal constants-- now these come from TkrGeometrySvc
    //double cal_top_z = -45.7;      // Meas off 1 Evt Disp.(29-may-03 - was -26.5) z co-ord. of top of Cal
    //double cal_half_width = 728.3; // measured from 1 Evt Disp - was (4*373.5(Tower Pitch) - 44(Cal Gap))/2
    //double cal_depth = 216.;       // Meas. off 1 Evt disp.. was 8 layers of the calorimeter

    // Cal constants, from detModel
    //double cal_top_z = m_calZTop;
    double cal_depth = -m_calZBot;
    double cal_half_width = 0.5*std::max(m_calXWidth, m_calYWidth);
    
    // Now do the leakage correction  
    // First: get the rad.lens. in the tracker 
    double t_tracker = track_1->getTkrCalRadlen();
    // Patch for error in KalFitTrack: 1/2 of first radiator left out
    int layer = track_1->getLayer();
    t_tracker += 0.5*m_geoSvc->getReconRadLenConv(layer)/costh;

    // Find the distance in Cal to nearest edge along shower axis
    //       Need to check sides as well as back
    Vector t_axis = axis.direction();     // This points in +z direction
    double s_top  = -(-m_calZTop + axis.position().z())/t_axis.z();
	double s_exit = -(axis.position().z())/t_axis.z(); 
    Point cal_top = axis.position(s_top); // Entry point into calorimeter
	Point tkr_exit= axis.position(s_exit);// Exit point from tracker

    CAL_z0    = cal_top.z();
    CAL_x0    = cal_top.x();
    CAL_y0    = cal_top.y();   

	// Only do leakage correction for tracks which "hit" the Calorimeter
	if(fabs(CAL_x0) > cal_half_width || fabs(CAL_y0) > cal_half_width) return sc; 


    double s_xp   = (-cal_half_width + tkr_exit.x())/t_axis.x();
    double s_xm   = ( cal_half_width + tkr_exit.x())/t_axis.x();
    double s_minx = (s_xp > s_xm) ? s_xp:s_xm; // Choose soln > 0. 
    
    double s_yp   = (-cal_half_width + tkr_exit.y())/t_axis.y();
    double s_ym   = ( cal_half_width + tkr_exit.y())/t_axis.y();
    double s_miny = (s_yp > s_ym) ? s_yp:s_ym; // Choose soln > 0. 
    
//    double s_minz = (m_calZTop - m_calZBot)/t_axis.z();
	double s_minz = (tkr_exit.z() - m_calZBot)/t_axis.z();
    // Now pick min. soln. of the x, y, and z sides 
    double s_min  = (s_minx < s_miny) ? s_minx:s_miny;  
    s_min         = (s_min  < s_minz) ? s_min :s_minz;
    
    // Set up a propagator to calc. rad. lens. 
	m_G4PropTool->setStepStart(tkr_exit, -t_axis);
	m_G4PropTool->step(s_min);  

	// Now loop over the steps to extract the materials
	int numSteps = m_G4PropTool->getNumberSteps();
	int istep  = 0;
	idents::VolumeIdentifier volId;
	idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();
	double radLen_CsI  = 0.;
	double radLen_Stuff = 0.;
	double radLen_Gap  = 0.; 
	double radLen_Cntr = 0.; 
	double radLen_CntrStuff = 0.; 
	double arcLen_CsI  = 0.; 
	double arcLen_Gap  = 0.; 
	double arcLen_Stuff = 0.; 
	double arcLen_Cntr = 0.; 
	for(; istep < numSteps; ++istep) {
		volId = m_G4PropTool->getStepVolumeId(istep);
		volId.prepend(prefix);
		double radLen_step = m_G4PropTool->getStepRadLength(istep);
		double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
		Point x_step       = m_G4PropTool->getStepPosition(istep);
		if(volId.size() > 8) { //This indicates being inside a CsI Xtal 
			radLen_CsI  += radLen_step;
			arcLen_CsI  += arcLen_step;
		}
		else {
			radLen_Stuff += radLen_step;
			arcLen_Stuff += arcLen_step;
			if(x_step.z() >= m_calZTop) {
				radLen_Gap += radLen_step;
				arcLen_Gap += arcLen_step;
			}
		}
		if(x_step.z() >= CAL_zEcntr) {
			radLen_Cntr += radLen_step;
			arcLen_Cntr += arcLen_step;
			if(volId.size() < 8) 
				radLen_CntrStuff += radLen_step;
		}
	}
    double t_cal_tot = radLen_CsI + radLen_Stuff; 
    double t_total   = t_tracker  + t_cal_tot;
    
    // Energy centroid in radiation lengths for this event.
    double t         = t_tracker + radLen_Cntr;

    // 2 Iterations to find energy
    // Note: the TMath incomplete gamma functions is really the fractional inner incomplete
	//       gamma function - i.e. its just what we need to do shower models.  
	//Leon: this is the algorithm that gave the best results for estimating
	//      leakage:
	//        a) In the shower model (wallet card) the ratio of the a & b paremters
	//           is the energy centroid (check it for yourself).  Also b is all but
	//           energy independent (not ours ... problem here -see below ... but the wallet cards)
	//           So given our event centroid and b - > make first guess for a.  But the ratio of 
	//           a/b = <t> only in an infinitely deep cal.  
	//        b) This gives first approx. for energy leakage fraction from model.
	//        c) Now find a & b from parametric  models parametric models using "corrected energy".
	//      I don't think is great. The reality is that first estimation correction raises
	//      the energy by ~ 10% at 10 GeV.  The second iteration using the a & b models gives
	//      a total correction of ~ + 33%. 
	//      The powerlaw parameterization of b and a are funny. Particularly b - our b for 
	//      CsI is ~ .46 and raising with E (i.e. too small - like smaller then uranium!  - and 
	//      changing too fast).  These are vistiges from the code in CalRecon. 
	//      Overall it seems to work albeit for heuristic methodology. Note however we may get strange
	//      results when used outside the energy range 50 MeV - 18 GeV.  I hope here that 
	//      the stuff Giebells is doing will come to the rescue - transistion from one to the other??
	//   
//  double a1 = alpha(ene_sum_corr);
    double b1 =  beta(ene_sum_corr); 
	double a1 = b1*t + 1.; //  + 1 is to compensates for finite cal effect.  This would better 
	                       // match it  to the final answer.  In reality it has a minor effect since all this 
	                       // typically affects things as log(E). 
    double in_frac_1 = TMath::Gamma(a1, b1*t_total);
    double a2 = alpha(ene_sum_corr/in_frac_1);
    double b2 =  beta(ene_sum_corr/in_frac_1); 
    double in_frac_2 = TMath::Gamma(a2, b2*t_total);
    
    ene_sum_corr /= in_frac_2; //Using 2nd order solution
    
    // Final Delta t factor (?).  Again - present in sumlated data - Needs to be understood 
	// Perhaps this is the event to event compensation for variations in the initiation of showers. 
    CAL_t_Pred      = a2/b2*(TMath::Gamma(a2+1.,b2*t_total)/
                             TMath::Gamma(a2,b2*t_total)); 
    CAL_deltaT      = t - CAL_t_Pred;
    
	// The "final" correction derived empirically from analyizing and flattening the 
	// resultant energy in cos(theta) and log10(E) 
	double logEsum = log10(ene_sum_corr); 
    //double ad_hoc_factor = (1.35-.15*logEsum )/(1.+.9*(.6 - costh)*(.6 - costh))/
		                       // (1.-.125*logEsum + (.125*logEsum)*costh);  
	//double ad_hoc_factor = (1.23 - .065*logEsum)/(1.+.14*(logEsum-1.)*(costh-.74));  
    double ad_hoc_factor = (1.08 + .02*logEsum);  
    // Store the results away 
    CAL_EneSum_Corr = ene_sum_corr * ad_hoc_factor;
    CAL_Energy_Corr = CAL_EnergySum*edge_corr * ad_hoc_factor/in_frac_2; 
    CAL_TotSum_Corr = CAL_EneSum_Corr/CAL_EnergySum;
    //CAL_Total_Corr  = CAL_Energy_Corr/CAL_EnergySum; 
	CAL_CsI_RLn     = radLen_CsI;
	CAL_MIP_Diff    = CAL_EnergySum - 12.07*radLen_CsI;
    CAL_Tot_RLn     = t_total;
    CAL_Cnt_RLn     = t;
    CAL_DeadTot_Rat = radLen_Stuff/t_total;
    CAL_DeadCnt_Rat = radLen_CntrStuff/t;
    //CAL_a_Parm      = a2;
    //CAL_b_Parm      = b2; 
    
    //CAL_Leak_Corr   = in_frac_1;
    CAL_Leak_Corr2  = in_frac_2;   
	CAL_TwrGap      = arcLen_Stuff - arcLen_Gap;

    return sc;
}

StatusCode CalValsTool::getCalInfo()
{
    m_calZTop = m_geoSvc->calZTop();
	m_calZBot = m_geoSvc->calZBot();
	m_calXWidth = m_geoSvc->calXWidth();
	m_calYWidth = m_geoSvc->calYWidth();

    return StatusCode::SUCCESS;
}
