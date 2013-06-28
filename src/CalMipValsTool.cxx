    
/** @file CalMipValsTool.cxx
@brief Calculates the Cal Mip analysis variables
@author F. Piron

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/CalMipValsTool.cxx,v 1.13 2013/06/26 01:05:29 lsrea Exp $
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
//@@@FP 07/09/05
#include "Event/Recon/CalRecon/CalMipClasses.h"
//@@@FP 07/09/05
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

/*! @class CalMipValsTool
@brief calculates Cal Mip values

  @author F. Piron
  */
  
class CalMipValsTool :   public ValBase
{
public:  
  CalMipValsTool( const std::string& type, 
                  const std::string& name, 
                  const IInterface* parent);
  
  virtual ~CalMipValsTool() { }
  
  StatusCode initialize();
  
  StatusCode calculate();
  
private:
  float m_num;
  float m_x0;
  float m_y0;
  float m_z0;

  double m_xDir;
  double m_yDir;
  double m_zDir;
  float m_d2edge;
  float m_arcLen;
  float m_ecor;
  float m_ecorRms;
  float m_chi2;
  float m_erm;
  //2007-04-05 Sylvain Guiriec
  float m_ermc;
  double m_derr;
  float m_barDist;
  //2007-04-05 Sylvain Guiriec End 
};
  
/** @page anatup_vars
    @section calmipvalstool CalMipValsTool Variables

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> CalMipNum  
<td>F<td>    Total number of found mip-like tracks in the Cal
<tr><td> CalMip[X/Y/Z]0  
<td>F<td>    [x/y/z] coordinates of the energy centroid of the best track
<tr><td> CalMip[X/Y/Z]Dir  
<td>D<td>    [x/y/z] direction cosines of the best track
<tr><td> CalMipD2edge  
<td>F<td>    Distance of the best track from the nearest edge of the Cal
<tr><td> CalMipArcLen  
<td>F<td>    Length of the best track (mm)
<tr><td> CalMipEcor  
<td>F<td>    Mean vertical-eq1uvalent energy (MeV) of the best track,
obtained by averaging the pathlength-corrected energies in each layer
<tr><td> CalMipEcorRms  
<td>F<td>    RMS of CalMipEcor
<tr><td> CalMipChi2  
<td>F<td>    Chi-squared of the direction fit for the best track
(combination of least squares in XZ and YZ planes)
<tr><td> CalMipErm 
<td>F<td>    total energy (MeV) contained in a cylinder of 1 Moliere radius 
around the best track
</table>

*/

// Static factory for instantiation of algtool objects
//static ToolFactory<CalMipValsTool> s_factory;
//const IToolFactory& CalMipValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(CalMipValsTool);
  
// Standard Constructor
CalMipValsTool::CalMipValsTool(const std::string& type, 
                               const std::string& name, 
                               const IInterface* parent)
               : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}
  
StatusCode CalMipValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    
    addItem("CalMipNum",     &m_num);      
    addItem("CalMipX0",      &m_x0);      
    addItem("CalMipY0",      &m_y0);      
    addItem("CalMipZ0",      &m_z0);      
    addItem("CalMipXDir",    &m_xDir);   
    addItem("CalMipYDir",    &m_yDir);    
    addItem("CalMipZDir",    &m_zDir);    
    addItem("CalMipD2edge",  &m_d2edge);  
    addItem("CalMipArcLen",  &m_arcLen);  
    addItem("CalMipEcor",    &m_ecor);    
    addItem("CalMipEcorRms", &m_ecorRms); 
    addItem("CalMipChi2",    &m_chi2);    
    addItem("CalMipErm",     &m_erm);

    //2007-04-05 Sylvain Guiriec
    addItem("CalMipErmc",    &m_ermc);
    addItem("CalMipDerr",    &m_derr);
    addItem("CalMipBarDist", &m_barDist);
    //2007-04-05 Sylvain Guiriec End 
    
    return sc;
}

StatusCode CalMipValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the calMipTrack collection   

    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << "FP : readCalMipTrackCol in CalMipValsTool" << endreq;

    SmartDataPtr<Event::CalMipTrackCol> p_calMipTrackCol(m_pEventSvc, EventModel::CalRecon::CalMipTrackCol); 

    Point  point(-1.,-1.,-1.);
    Vector dir(-1.,-1.,-1.);
    double d2c=-6.;
    double d2edge=-6.;
    int    calEdge=-6;
    double arcLen=-6.;
    double ecor=-6.;
    double ecorRms=-6.;
    double chi2=-6.;
    double erm=-6.;
    //100 tracks maxi, and 2 criteria :
    //   1/ at the CAL boundary
    //   2/ Ecor in energy interval given by 25% of MIP Landau function maximum
    int ipassed[100][2];
    for(int j=0;j<100;j++)
    {
        ipassed[j][0]=0;
        ipassed[j][1]=0;
    }
    int n1=0,n2=0;

    m_num=0;
    m_x0=-5.;
    m_y0=-5.;
    m_z0=-5.;
    m_xDir=-5.;
    m_yDir=-5.;
    m_zDir=-5.;
    m_d2edge=-5.;
    m_arcLen=-5.;
    m_ecor=-5.;
    m_ecorRms=-5.;
    m_chi2=-5.;
    m_erm=-5.;

    if (p_calMipTrackCol==0)
      return sc;
    m_num=p_calMipTrackCol->size();
    int it=-1;
    for(Event::CalMipTrackCol::const_iterator calMipTrackIter=p_calMipTrackCol->begin(); calMipTrackIter != p_calMipTrackCol->end(); calMipTrackIter++)
    {
        it++;
        log << "------------------------------------------------------------" << endreq;
        log << "counter col=" << it << " / size = " << m_num << endreq;
        Event::CalMipTrack* p_calMipTrack    =  *calMipTrackIter;
        p_calMipTrack->writeOut(log);
        log << "------------------------------------------------------------" << endreq;

        point   = p_calMipTrack->getPoint   ();
        dir     = p_calMipTrack->getDir     ();
        d2c     = p_calMipTrack->getD2C     ();
        d2edge  = p_calMipTrack->getD2Edge  ();
        calEdge = p_calMipTrack->getCalEdge ();
        arcLen  = p_calMipTrack->getArcLen  ();
        ecor    = p_calMipTrack->getEcor    ();
        ecorRms = p_calMipTrack->getEcorRms ();
        chi2    = p_calMipTrack->getChi2    ();
        erm     = p_calMipTrack->getErm     ();
            
        if (d2edge<1000000)
        {
            n1++;
            ipassed[it][0]=1;
            if (ecor>10.0 && ecor<14.2)
            {
                n2++;
                ipassed[it][1]=1;
            }
        }
    }
        
//     if (n2<=0)
    if (n1<=0)
      return sc;
    else
    {
        // loop again over tracks to keep the track with best chi2
        int it=-1;
        double chi2min=999999;
        for(Event::CalMipTrackCol::const_iterator calMipTrackIter=p_calMipTrackCol->begin(); calMipTrackIter != p_calMipTrackCol->end(); calMipTrackIter++)
        {
            it++;
            //        if (!ipassed[it][1])
            if (!ipassed[it][0])
                continue;
            
            Event::CalMipTrack* p_calMipTrack    =  *calMipTrackIter;
        
            double chi2Cur    = p_calMipTrack->getChi2    ();
            if (chi2Cur<chi2min)
            {
                chi2min=chi2Cur;
                chi2=chi2Cur;
                point   = p_calMipTrack->getPoint   ();
                dir     = p_calMipTrack->getDir     ();
                d2c     = p_calMipTrack->getD2C     ();
                d2edge  = p_calMipTrack->getD2Edge  ();
                calEdge = p_calMipTrack->getCalEdge ();
                arcLen  = p_calMipTrack->getArcLen  ();
                ecor    = p_calMipTrack->getEcor    ();
                ecorRms = p_calMipTrack->getEcorRms ();
                erm     = p_calMipTrack->getErm     ();
            }
        }
    }
    
    m_x0=point.x();
    m_y0=point.y();
    m_z0=point.z();
    m_xDir=dir.x();
    m_yDir=dir.y();
    m_zDir=dir.z();
    m_d2edge=d2edge;
    m_arcLen=arcLen;
    m_ecor=ecor;
    m_ecorRms=ecorRms;
    m_chi2=chi2;
    m_erm=erm;

    //2007-04-05 Sylvain Guiriec

    double ermc=0.;
    double derr=0.;
    double barDist=0.;

    m_ermc=-100.;
    m_derr=-100.;
    m_barDist=-100.;

    //Recover EventHeader Pointer
    //SmartDataPtr<Event::EventHeader> pEvent(m_pEventSvc, EventModel::EventHeader);

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrTrackCol>   pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol);
    if(!pTracks) return sc;

    //SmartDataPtr<Event::TkrVertexCol>  pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    //SmartDataPtr<Event::TkrClusterCol> pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    // all variable values are preset to zero. Be sure to re-initialize the ones you care about  

    //double die_width = m_tkrGeom->ladderPitch();
    //int nDies = m_tkrGeom->nWaferAcross();

    if (!pTracks || n1<=0)
      return sc;
    else
      {
        // Count number of tracks
        int nTracks = pTracks->size();
        //int Tkr_No_Tracks   = nTracks;
        
        if(nTracks < 1) return sc;
        
        // Get the first Track - it should be the "Best Track"
        Event::TkrTrackColConPtr pTrack = pTracks->begin();
        
        const Event::TkrTrack* track_1 = *pTrack;
        
        Point  x1 = track_1->getInitialPosition();
        Vector t1 = track_1->getInitialDirection();
        
        double Tkr_1_xdir        = t1.x();
        double Tkr_1_ydir        = t1.y();
        double Tkr_1_zdir        = t1.z();
        
        float Tkr_1_x0          = x1.x();
        float Tkr_1_y0          = x1.y();
        float Tkr_1_z0          = x1.z();

        if (m_arcLen>0)
          ermc=m_erm/m_arcLen;
        else if (m_arcLen<=0)
          ermc=-100.;
        //@@@@@
        derr=TMath::ACos(m_xDir*Tkr_1_xdir+
                                m_yDir*Tkr_1_ydir+
                                m_zDir*Tkr_1_zdir);
        // could just be:
        //derr = TMath::ACos(dir.dot(t1));

        
        //@@@@@ Distance between tracker track & calorimeter track's centroid point
        double l=(Tkr_1_xdir*(m_x0-Tkr_1_x0)+Tkr_1_ydir*(m_y0-Tkr_1_y0)+Tkr_1_zdir*(m_z0-Tkr_1_z0))/(TMath::Power(Tkr_1_xdir,2)+TMath::Power(Tkr_1_ydir,2)+TMath::Power(Tkr_1_zdir,2));
        double Bx   = l*Tkr_1_xdir+Tkr_1_x0;
        double By   = l*Tkr_1_ydir+Tkr_1_y0;
        double Bz   = l*Tkr_1_zdir+Tkr_1_z0;
        barDist = TMath::Sqrt(TMath::Power(m_x0-Bx,2)+TMath::Power(m_y0-By,2)+TMath::Power(m_z0-Bz,2));
      }

    m_ermc=ermc;
    m_derr=derr;
    m_barDist=barDist;

    //2007-04-05 Sylvain Guiriec End  

    log << MSG::DEBUG << " FP : readCalMipTrackCol in CalMipValsTool End  " << m_num << endreq;
    return sc;
}

