/** @file CalValsTool.cxx
@brief Calculates the Cal analysis variables
@author Bill Atwood, Leon Rochester

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/CalValsTool.cxx,v 1.15 2003/05/08 15:51:30 atwood Exp $
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
//#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TMath.h"

#ifndef M_PI
#define M_PI = 3.14159265358979323846
#endif

namespace {
    const int _nTowers = 4;
    
    double signC(double x) { return (x<0.) ? -1.:1.;} 
    
    double twrEdgeC(double x, double y, double pitch, int &XY, int &outer) {
        double edge = 0.; 
        double x_twr = signC(x)*(fmod(fabs(x),pitch) - pitch/2.);
        double y_twr = signC(y)*(fmod(fabs(y),pitch) - pitch/2.);
        
        outer = 0; 
        
        if(fabs(x_twr) > fabs(y_twr)) {
            edge = pitch/2. - fabs(x_twr);
            XY = 1; 
            if(fabs(x) > 0.5*(_nTowers-1)*pitch) outer = 1;
        }
        else {
            edge = pitch/2. - fabs(y_twr);
            XY = 2;
            if(fabs(y) > 0.5*(_nTowers-1)*pitch) outer = 1;
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
    
    double contained_frac(double x, double y, double pitch, double gap,  
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
        if(fabs(x) < 1.5*pitch) { // X edge is not outside limit of LAT
            double r_frac_minus = (edge + gap_x/2.)/r;
            in_frac_x += circle_frac_simp(-r_frac_minus, angle_factor);
        }
        
        // Y Edges
        double y_twr = signC(y)*(fmod(fabs(y),pitch) - pitch/2.);
        edge = pitch/2. - fabs(y_twr);
        r_frac_plus = (edge-gap_y/2.)/r; 
        angle_factor = cos(phi)*(1./costh - 1.);
        double in_frac_y  =  circle_frac_simp(r_frac_plus, angle_factor);
        if(fabs(y) < 1.5*pitch) { // X edge is not outside limit of LAT
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
        
        // double b      = exp(0.031*log(E/1000.))/1.95;  //RR175
        double b      = exp(0.031*log(E/1000.))/2.3;    //RR176
        return b;
    }
    
    double alpha(double E) {
        // Following expression for a (alpha) comes from CalRecon
        // double alpha      = 2.65*exp(0.15*log(E/1000.));
        // The above expression gives too large values
        
        //  double a = 2.75*exp(0.14*log(E/1000.));  //RR175
        //  double a = 2.90*exp(0.115*log(E/1000.)); //RR176
        //  double a = 2.70*exp(0.15*log(E/1000.));  //RR179
        //  double a = 2.50*exp(0.18*log(E/1000.));  //RR180
        double a = 2.50*exp(0.19*log(E/1000.));    //RR182
        return a;
    }
    
    double gamma_alpha(double a, double b) {
        double val = TMath::Gamma(a+1., b)/TMath::Gamma(a, b); 
        return val;
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
      
      IKalmanParticle*       pKalParticle;   
      /// the GlastDetSvc used for access to detector info
      IGlastDetSvc*    m_detSvc; 
      /// store the towerPitch
      double m_towerPitch;
      /// 
      IPropagatorSvc* m_propSvc;
      
      //Global Calorimeter Tuple Items
      double CAL_EnergySum;
      double CAL_Leak_Corr; 
      double CAL_Leak_Corr2;
      double CAL_Edge_Corr; 
      double CAL_EdgeSum_Corr;     
      double CAL_Total_Corr; 
      double CAL_TotSum_Corr; 
      double CAL_Tot_RLn;
      double CAL_Cnt_RLn; 
      double CAL_DeadTot_Rat;
      double CAL_DeadCnt_Rat; 
      double CAL_a_Parm;
      double CAL_b_Parm; 
      double CAL_t_Pred; 
      double CAL_deltaT;
      
      double CAL_EneSum_Corr;
      double CAL_Energy_Corr;
      double CAL_x0;
      double CAL_y0;
      double CAL_z0;
      double CAL_xdir;
      double CAL_ydir;
      double CAL_zdir;
      double CAL_TwrEdge;
      double CAL_TE_Nrm;
      double CAL_Track_Sep;
	  double CAL_Xtal_Ratio;
	  double CAL_eLayer[8];
	  double CAL_No_Xtals_Trunc;
	  double CAL_Long_Rms;
	  double CAL_Trans_Rms;

      
      //Calimeter items with Recon - Tracks
      double CAL_Track_DOCA;
      double CAL_Track_Angle;
	  double CAL_TwrGap; 
      double CAL_x0_corr;
      double CAL_y0_corr;
      double CAL_z0_corr;
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
          
          // find the GlastDevSvc service
          if (service("GlastDetSvc", m_detSvc).isFailure()){
              log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
              return StatusCode::FAILURE;
          }
          m_detSvc->getNumericConstByName("towerPitch", &m_towerPitch);
                   
          // pick up the chosen propagator
          if (service("GlastPropagatorSvc", m_propSvc).isFailure()) {
              log << MSG::ERROR << "Couldn't find the GlastPropagatorSvc!" << endreq;
              return StatusCode::FAILURE;
          }
          
          pKalParticle = m_propSvc->getPropagator();
          
      } else {
          return StatusCode::FAILURE;
      }
      
      // load up the map
      
      addItem("CalEnergySum",  &CAL_EnergySum);
      addItem("CalEnergyCorr", &CAL_Energy_Corr);
      addItem("CalEneSumCorr", &CAL_EneSum_Corr);
      
      addItem("CalLeakCorr",   &CAL_Leak_Corr);
      addItem("CalLeakCorr2",  &CAL_Leak_Corr2);
      
      addItem("CalEdgeCorr",   &CAL_Edge_Corr);
      addItem("CalEdgeSumCorr",&CAL_EdgeSum_Corr);
      addItem("CalTotalCorr",  &CAL_Total_Corr);
      addItem("CalTotSumCorr", &CAL_TotSum_Corr);
      
      addItem("CalTotRLn",     &CAL_Tot_RLn);
      addItem("CalCntRLn",     &CAL_Cnt_RLn);
      addItem("CalDeadTotRat", &CAL_DeadTot_Rat);
      addItem("CalDeadCntRat", &CAL_DeadCnt_Rat);
      addItem("CalAParm",      &CAL_a_Parm);
      addItem("CalBParm",      &CAL_b_Parm);
      addItem("CalTPred",      &CAL_t_Pred);
      addItem("CalDeltaT",     &CAL_deltaT);
      
      addItem("CalTwrEdge",    &CAL_TwrEdge);
      addItem("CalTENrm",      &CAL_TE_Nrm);
      addItem("CalTrackSep",   &CAL_Track_Sep);
      addItem("CalTrackDoca",  &CAL_Track_DOCA);
	  addItem("CalTwrGap",     &CAL_TwrGap);
      addItem("CalXtalRatio",  &CAL_Xtal_Ratio);
	  addItem("CalELayer0",    &CAL_eLayer[0]);
      addItem("CalELayer1",    &CAL_eLayer[1]);
	  addItem("CalELayer2",    &CAL_eLayer[2]);
	  addItem("CalELayer3",    &CAL_eLayer[3]);
	  addItem("CalELayer4",    &CAL_eLayer[4]);
	  addItem("CalELayer5",    &CAL_eLayer[5]);
	  addItem("CalELayer6",    &CAL_eLayer[6]);
	  addItem("CalELayer7",    &CAL_eLayer[7]);
	  addItem("CalXtalsTrunc", &CAL_No_Xtals_Trunc);
	  addItem("CalLongRms",    &CAL_Long_Rms);
	  addItem("CalTransRms",   &CAL_Trans_Rms);

      addItem("CalX0",         &CAL_x0);
      addItem("CalY0",         &CAL_y0);
      addItem("CalZ0",         &CAL_z0);
      addItem("CalXDir",       &CAL_xdir);
      addItem("CalYDir",       &CAL_ydir);
      addItem("CalZDir",       &CAL_zdir);
      
      addItem("CalX0Corr",     &CAL_x0_corr);
      addItem("CalY0Corr",     &CAL_y0_corr);
      addItem("CalZ0Corr",     &CAL_z0_corr);

      
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

	CAL_Long_Rms = calCluster->getRmsLong();
	CAL_Trans_Rms = calCluster->getRmsTrans();

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
    
    CAL_x0          = cal_pos.x();
    CAL_y0          = cal_pos.y();
    CAL_z0          = cal_pos.z();
    CAL_xdir        = cal_dir.x();
    CAL_ydir        = cal_dir.y();
    CAL_zdir        = cal_dir.z();
    
    double twr_pitch = m_towerPitch;
    int iView = 0;
    int outside = 0;
    CAL_TwrEdge = twrEdgeC(CAL_x0, CAL_y0, twr_pitch, iView, outside);
    
    if(iView==1) CAL_TE_Nrm  = cal_dir.x();
    else          CAL_TE_Nrm  = cal_dir.y();
    
    // with no tracks we've done what we can do.
    if(!pTracks) return sc; 
    int num_tracks = pTracks->size(); 
    if(num_tracks <= 0 ) return sc;
    
    // Get the first track
    Event::TkrFitConPtr pTrack1 = pTracks->begin();
    Event::TkrKalFitTrack* track_1  
        = dynamic_cast<Event::TkrKalFitTrack*>(*pTrack1);
    
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
    Ray axis(cal_pos, x_diff.unit());
    double arc_len = cal_pos.z()/t0.z(); 
    Point z_0 = axis.position(arc_len); 
    axis = Ray(z_0, x_diff.unit());
    
    std::vector<Vector> pos_Layer = calCluster->getPosLayer();
    std::vector<double> ene_Layer = calCluster->getEneLayer();
    
    double ene_sum_corr = 0.;
    Vector pos_corr(0., 0., 0.); 
    
    // Factors controlling edge correction: core radius, fringe radius, core fraction
    //    and size of inter-tower gap (active-to-active area) 
    
    double rm_hard   = 18.; // Fixed at ~ 1 rad. len in CsI 
    
    //  double rm_soft   = 55.+ 200.*cal_trans((CAL_EnergySum-120)/110.); 
    double rm_soft   = 55.+ 200.*cal_trans((CAL_EnergySum-250)/200.); 
    
    // double hard_frac =.50 +.35*cal_trans((670.- CAL_EnergySum)/200.); 
    double hard_frac =.50 +.35*cal_trans((670.- CAL_EnergySum)/300.);   
    
    
    double gap       = 48;// + 30.*cal_trans((CAL_EnergySum-80)/70.);
    
    // Reset shower area circle when multiple tracks are present
    // This increasingly matters below 1000 MeV.
    if(num_tracks > 1) {
        pTrack1++;
        Event::TkrKalFitTrack* track_2  
            = dynamic_cast<Event::TkrKalFitTrack*>(*pTrack1);
        Event::TkrFitTrackBase::TrackEnd end = Event::TkrFitTrackBase::End; 
        Ray trj_1 = Ray(track_1->getPosition(end), track_1->getDirection(end)); 
        Ray trj_2 = Ray(track_2->getPosition(end), track_2->getDirection(end));
        double delta_z = trj_1.position().z() - (-26.5);
        double arc_len = delta_z/fabs(trj_1.direction().z());
        Point cal_1 = trj_1.position(arc_len); 
        delta_z = trj_2.position().z() - (-26.5);
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
            twr_pitch, gap, rm_soft, costh, phi_90);
        
        // Cut off correction upon leaving (through a side)
        if(in_frac_soft < .5) break;
        
        double in_frac_hard = contained_frac(xyz_layer.x(), xyz_layer.y(), 
            twr_pitch, gap, rm_hard, costh, phi_90);
        
        double corr_factor 
            = 1./((1.-hard_frac)*in_frac_soft + hard_frac*in_frac_hard);
        if(corr_factor > max_corr) corr_factor = max_corr;  
        double ene_corr = ene_Layer[i]*corr_factor;
        ene_sum_corr += ene_corr;
        pos_corr     += ene_corr*pos_Layer[i]; 
        edge_corr    += corr_factor;
        good_layers  += 1.; 
    }
    if(ene_sum_corr < 1.) return sc;
    
    pos_corr /= ene_sum_corr;
    CAL_z0_corr = pos_corr.z();
    CAL_x0_corr = pos_corr.x();
    CAL_y0_corr = pos_corr.y();
    CAL_EdgeSum_Corr = ene_sum_corr/CAL_EnergySum;
    if (good_layers>0) edge_corr /= good_layers;
    CAL_Edge_Corr = edge_corr; 
    
    // Set some Calorimeter constants - should come from detModel
    double cal_top_z = -26.5;      // z co-ord. of top of Cal
    double x0_CsI    = 18.5;       // Rad. Len. of CsI
    double x0_Crap   = 100.;       // Rad Len. of junk inactive junk in Cal
    double xtal_dz   = 19.9;       // Xtal height in z
    double cal_half_width = 725.;  // (4*373.5(Tower Pitch) - 44(Cal Gap))/2
    double cal_depth = 8*xtal_dz;  // 8 layers of in the calorimeter
    
    // Now do the leakage correction  
    // First: get the rad.lens. in the tracker 
    double t_tracker = track_1->getTkrCalRadlen();
    // Patch for error in KalFitTrack: 1/2 of first radiator left out
    if(track_1->getLayer() < 12) t_tracker += .015/costh; 
    else                         t_tracker += .09/costh; 
    
    // Find the distance in Cal to nearest edge along shower axis
    // Need to check sides as well as back
    Vector t_axis = axis.direction();     // This points in +z direction
    double s_top  = -(-cal_top_z + axis.position().z())/t_axis.z(); 
    Point cal_top = axis.position(s_top); // Entry point into calorimeter
    
    double s_xp   = (-cal_half_width + cal_top.x())/t_axis.x();
    double s_xm   = ( cal_half_width + cal_top.x())/t_axis.x();
    double s_minx = (s_xp > s_xm) ? s_xp:s_xm; // Choose soln > 0. 
    
    double s_yp   = (-cal_half_width + cal_top.y())/t_axis.y();
    double s_ym   = ( cal_half_width + cal_top.y())/t_axis.y();
    double s_miny = (s_yp > s_ym) ? s_yp:s_ym; // Choose soln > 0. 
    
    double s_minz = cal_depth/t_axis.z();
    // Now pick min. soln. of the x, y, and z sides 
    double s_min  = (s_minx < s_miny) ? s_minx:s_miny;  
    s_min         = (s_min  < s_minz) ? s_min :s_minz;
    
    // Set up a propagator to calc. rad. lens. 
    pKalParticle->setStepStart(cal_top, -t_axis, s_min);
    
    double t_cal_tot  = pKalParticle->radLength();
    double t_cal_tot2 = s_min/(1.063*x0_CsI);
    double t_total    = t_tracker + t_cal_tot;
    double t_total2   = t_tracker + t_cal_tot2; 
    double t_cal_dead = (s_min - t_cal_tot*x0_CsI)/(1.-x0_CsI/x0_Crap)
        /x0_Crap;
    
    // Energy centroid in radiation lengths for this event in rad.lens.
    double s_t_cal   = (cal_top_z-CAL_z0)/costh; 
    double t_cal_cnt = pKalParticle->radLength(s_t_cal);
    double t_cal_cnt2= s_t_cal/(1.063*x0_CsI); 
    double t         = t_tracker + t_cal_cnt; 
    double t2        = t_tracker + t_cal_cnt2; 
    double t_dead    = (s_t_cal - t_cal_cnt*x0_CsI)/(1.-x0_CsI/x0_Crap)
        /x0_Crap; 
    
    // Following expression for b (beta) comes from CalRecon
    // double beta      = exp(0.031*log(CAL_EnergySum/1000.))/2.29;
    // The above expression gives too small values
    //double beta0      = exp(0.061*log(ene_sum_corr/1000.))/2.25; //RR177
    double beta0      = exp(0.081*log(ene_sum_corr/1000.))/2.25; //RR178
    
    double beta_t  = beta0 * t_total2;
    double alpha_t = beta0 * t2;
    if(t > 1.2 * t_total) {
        alpha_t = 1.2*beta_t; 
    } 
    else {
        // Iterations solution for alpha
        alpha_t   /= gamma_alpha(alpha_t, beta_t);
        double alpha_t2   = alpha_t/(gamma_alpha(alpha_t, beta_t));
        alpha_t  = .8*alpha_t + .2*alpha_t2; 
    }
    // Don't accept giant leakage corrections
    double in_frac = TMath::Gamma(alpha_t, beta_t); 
    if(in_frac < .2) {  
        in_frac = .2; 
    }
    // 2 Iterations to find energy
    double a1 = alpha(ene_sum_corr);
    double b1 = beta(ene_sum_corr); 
    double in_frac_1 = TMath::Gamma(a1, b1*t_total);
    double a2 = alpha(ene_sum_corr/in_frac_1);
    double b2 =  beta(ene_sum_corr/in_frac_1); 
    double in_frac_2 = TMath::Gamma(a2, b2*t_total);
    
    ene_sum_corr /= in_frac; //Using <t> solution 
    
    // Final cos(theta) dependence (?)
    double cos_factor = .75 + .25*costh;
    
    // Final Delta t factor (?)
    CAL_t_Pred      = a2/b2*(TMath::Gamma(a2+1.,b2*t_total)/
        TMath::Gamma(a2,b2*t_total)); 
    double delta_t = t - CAL_t_Pred; 
    double dt_factor  = 1.04/(1.+.04*delta_t); // 1 GeV -> .0073
    
    ene_sum_corr *= dt_factor/cos_factor;
    
    // Store the results away 
    CAL_EneSum_Corr = ene_sum_corr;
    CAL_Energy_Corr = CAL_EnergySum*edge_corr*dt_factor/in_frac/cos_factor; 
    CAL_TotSum_Corr = CAL_EneSum_Corr/CAL_EnergySum;
    CAL_Total_Corr  = edge_corr/in_frac/cos_factor; 
    CAL_Tot_RLn     = t_total;
    CAL_Cnt_RLn     = t;
    CAL_DeadTot_Rat = t_cal_dead/t_total;
    CAL_DeadCnt_Rat = t_dead/t;
    CAL_a_Parm      = a2;
    CAL_b_Parm      = b2; 
    
    CAL_Leak_Corr   = in_frac;
    CAL_Leak_Corr2  = in_frac_2;   
    CAL_deltaT      = CAL_Cnt_RLn - CAL_t_Pred;

	    // New section to compute arclength of trajectory (caltop - taxis) falling between towers
	double x_twr = signC(-t_axis.x())*(fmod(fabs(cal_top.x()),twr_pitch) - twr_pitch/2.);
    double y_twr = signC(-t_axis.y())*(fmod(fabs(cal_top.y()),twr_pitch) - twr_pitch/2.);
    double x_slope   = (fabs(t_axis.x()) > .0001)? -t_axis.x():.00001;
    double s_x       = (signC(x_slope)*twr_pitch/2. - x_twr)/x_slope; 
	double y_slope   = (fabs(t_axis.y()) > .0001)? -t_axis.y():.00001;
    double s_y       = (signC(y_slope)*twr_pitch/2. - y_twr)/y_slope;
		
	CAL_TwrGap = 0.; 
    if(s_x < s_y) { // Goes out x side of CAL Module
		if(cal_top.z() - s_x*t_axis.z() > -160) {
			double s_max = (160.+cal_top.z())/t_axis.z();
			CAL_TwrGap = gap/fabs(x_slope);
			if((CAL_TwrGap + s_x)> s_max ) CAL_TwrGap = s_max-s_x;
		}
	}
	else {          // Goes out y side
		if(cal_top.z() - s_y*t_axis.z() > -160) {
			double s_max = (160.+cal_top.z())/t_axis.z();
			CAL_TwrGap = gap/fabs(y_slope);
			if((CAL_TwrGap + s_y)> s_max ) CAL_TwrGap = s_max-s_y;
		}
	}
    
    return sc;
}
