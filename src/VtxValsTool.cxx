
/** @file VtxValsTool.cxx
@brief Calculates the Vtx analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/VtxValsTool.cxx,v 1.30.366.1 2010/09/18 02:35:15 heather Exp $
*/

// Include files


#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "Doca.h"

namespace {
	double erf(double x) {
		double t = 1./(1.+.47047*x);
		double results = 1. -(.34802 - (.09587 -.74785*t)*t)*t*exp(-x*x);
		return results;
	}
	double thrshold(double x) {
		if(x < 0) return (.5*(1. + erf(-x)));
		else      return (.5*(1. - erf( x)));
	}

}

/*! @class VtxValsTool
@brief calculates Vtx values

@authors Bill Atwood, Leon Rochester
*/

class VtxValsTool : public ValBase
{
public:

	VtxValsTool( const std::string& type, 
		const std::string& name, 
		const IInterface* parent);

	virtual ~VtxValsTool() { }

	StatusCode initialize();

	StatusCode calculate();

private:

	//Vertexing Items
    int   VTX_numVertices;
	float VTX_xdir;
	float VTX_ydir;
	float VTX_zdir;
	float VTX_Phi;
	float VTX_Theta;
    float VTX_Sxx;
    float VTX_Sxy;
    float VTX_Syy;
    float VTX_ThetaErr;
    float VTX_PhiErr;
    //float VTX_ErrAsym;
    //float VTX_CovDet;
	float VTX_x0;
	float VTX_y0;
	float VTX_z0;
	float VTX_Angle;
	float VTX_DOCA;
	float VTX_Head_Sep;

	float VTX_S1;
	float VTX_S2;

	float VTX_Quality; 
	float VTX_Chisq; 
	float VTX_AddedRL;
	float VTX_Status;

    float VTX2_transDoca;
    float VTX2_longDoca;

    float VTX2_xdir;
	float VTX2_ydir;
	float VTX2_zdir;
	float VTX2_Phi;
	float VTX2_Theta;
	float VTX2_x0;
	float VTX2_y0;
	float VTX2_z0;
	float VTX2_Angle;
	float VTX2_DOCA;
	float VTX2_Head_Sep;
    float VTX2_Status;

    float VTXN_xdir;
	float VTXN_ydir;
	float VTXN_zdir;

	float VTXN_Sxx;
    float VTXN_Sxy;
    float VTXN_Syy;
	float VTXN_ChgWt;
	float VTXN_NeutWt;

	float VTXN1_xdir;
	float VTXN1_ydir;
	float VTXN1_zdir;

	float VTXN1_Sxx;
    float VTXN1_Sxy;
    float VTXN1_Syy;
    float VTXN1_ChgWt;
	float VTXN1_NeutWt;
};

// Static factory for instantiation of algtool objects
//static ToolFactory<VtxValsTool> s_factory;
//const IToolFactory& VtxValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(VtxValsTool);

// Standard Constructor
VtxValsTool::VtxValsTool(const std::string& type, 
						 const std::string& name, 
						 const IInterface* parent)
						 : ValBase( type, name, parent )
{    
	// Declare additional interface
	declareInterface<IValsTool>(this); 
}

/** @page anatup_vars 
    @section vtxvalstool VtxValsTool Variables

In what follows below, whenever the first 2 vertices are referenced, they will be called "Vtx" and
"Vtx2". (This is denoted by "Vtx[/2]Xxx".) This is to maintain backward compatibility with existing code.
If there is no second vertex in the event, the Vtx2 quanitities are set to zero.

We've added some variables associated with a new kind of vertex: the "neutral" vertex.
This vertex is made by combining the direction of a simple "charged" vertex with the direction of the line
between the vertex point and the centroid of the CAL energy, weighted by the covariance matrices and by
an empirical factor involving the energy, chi-squareds, opening angle, etc. In the table below, "VtxNeut"
refers to the neutral vertex made from the first ("best") charged vertex, and "VtxNeut1" to a the neutral vertex made
from the best track. (For 1-track vertices, these would be the same.)

("Vtx[/Neut/Neut1]XXX" means that there are 3 variables: VtxXxx, VtxNeutXxx, and VtxNeut1Xxx.)

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> VtxNumVertices
<td>I<td>   Number of vertices in the event
<tr><td> Vtx[/2/Neut/Neut1][X/Y/Z]Dir 
<td>F<td>   [x/y/z] direction cosine of the vertex: 
            the first is "Vtx"; the 2nd is "Vtx2".
            For Neut and Neut1, see above. 
<tr><td> VtxPhi 
<td>F<td>   Azimuthal angle of vertex, radians 
            (direction of source, not flight direction!) Range: (0,2pi) 
<tr><td> VtxTheta 
<td>F<td>   Polar angle of vertex, radians (ditto direction)
<tr><td> VtxThetaErr  
<td>F<td>   Error on the measurement of theta  
<tr><td> VtxPhiErr  
<td>F<td>   Error on the measurement of phi.  
<tr><td> Vtx[/Neut/Neut1]S[XX/YY]  
<td>F<td>   [x-x/y-y] element of the covariance matrix; square of error on [x/y]  
<tr><td> Vtx[/Neut/Neut1]SXY  
<td>F<td>   x-y element of the covariance matrix; covariance  
<tr><td> Vtx[/2][X/Y/Z]0 
<td>F<td>   [x/y/z] coordinate of vertex The first vertex is "Vtx"; 
            the 2nd is "Vtx2". If the two tracks making up 
            the vertex are nearly parallel, 
            the coordinates of the vertex may become very large. 
<tr><td> Vtx2TransDoca
<tr><td> Vtx[/2]Angle 
<td>F<td>   Angle between the two tracks of the vertex (radians);             
            the first vertex is "Vtx"; the 2nd is "Vtx2".  
<tr><td> Vtx[/2]DOCA 
<td>F<td>   Distance of closest approach between the two tracks; 
            the first vertex is "Vtx"; the 2nd is "Vtx2". 
<tr><td> Vtx[/2]HeadSep 
<td>F<td>   Distance between the heads of the two tracks; 
            the first vertex is "Vtx"; the 2nd is "Vtx2".
<tr><td> Vtx2LongDoca
<td>F<td>   Longitudinal distance of the 2nd vertex along the axis of the first vertex;
            positive if the 2nd vertex is below the first.
<tr><td> Vtx2TransDoca
<td>F<td>   Transverse distance of the 2nd vertex point from the axis of the first vertex
<tr><Td> Vtx[/2]Status
<td>F<td>   Summary of track composition and topology.
            See TkrVertex.h in the Event package for the current description.
            The definitions as of GR v7r2 are:
@verbatim

    |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
                                         [Track Topology ] [Track composition]

    ONETKRVTX = 0x0001,  //Set if single track vertex
    TWOTKRVTX = 0x0002,  //Set if 2 track vertex
    MULTKRVTX = 0x0004,  //Set if >2 track vertex
	NEUTRALVTX = 0x0008,  //Set if vertex includes neutral energy vector
    DOCAVTX   = 0x0010,  //Set if vertex location set by DOCA point
    FIRSTHIT  = 0x0020,  //Set if two tracks share first hit
    STAGVTX   = 0x0040,  //Set if tracks don't start in same plane (staggered)
    CROSSTKR  = 0x0080   //Set if DOCA location lies inside track hits

@endverbatim
<tr><td> VtxQuality 
<td>F<td>   Vertex quality parameter used to order the possible vertices and 
            select the best one. <strong>Should generally not be used in analysis</strong>. 
<tr><td> VtxChisq 
<td>F<td>    The covariant chi-squared for the pairing of the tracks. 
<tr><td> VtxS[1/2] 
<td>F<td>   Distance of DOCA point from head of track [1/2] 
<tr><td> VtxAddedRL 
<td>F<td>   The additional radiation lengths prior to the first measured silicon 
            strip hit at the vertex location. <em>New!</em> 
</table>
*/


StatusCode VtxValsTool::initialize()
{
	StatusCode sc = StatusCode::SUCCESS;

	MsgStream log(msgSvc(), name());

	if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

	// get the services    

	if( serviceLocator() ) {        
	} else {
		return StatusCode::FAILURE;
	}

	// load up the map

	// Pair reconstruction
    addItem("VtxNumVertices", &VTX_numVertices);
	addItem("VtxXDir",      &VTX_xdir);     
	addItem("VtxYDir",      &VTX_ydir);     
	addItem("VtxZDir",      &VTX_zdir);     
	addItem("VtxPhi",       &VTX_Phi);  
	addItem("VtxTheta",     &VTX_Theta);  

    addItem("VtxThetaErr",   &VTX_ThetaErr);
    addItem("VtxPhiErr",     &VTX_PhiErr);
    addItem("VtxSXX",        &VTX_Sxx);
    addItem("VtxSXY",        &VTX_Sxy);
    addItem("VtxSYY",        &VTX_Syy);

    /* in case we want this later
    <tr><td> VtxErrAsym  
    <td>F<td>   Tkr1SXY/(Tkr1SXX + Tkr1SYY)  
    <tr><td> VtxCovDet  
    <td>F<td>   Determinant of the error matrix, 
             but normalized to remove the dependence on cos(theta)
    */
    //addItem("VtxErrAsym",    &VTX_ErrAsym);
    //addItem("VtxCovDet",     &VTX_CovDet);

   
    addItem("VtxX0",        &VTX_x0);       
	addItem("VtxY0",        &VTX_y0);       
	addItem("VtxZ0",        &VTX_z0);

	addItem("VtxAngle",     &VTX_Angle);    
	addItem("VtxDOCA",      &VTX_DOCA);
	addItem("VtxHeadSep",   &VTX_Head_Sep);
	addItem("VtxStatus",    &VTX_Status); 
	addItem("VtxQuality",   &VTX_Quality);
	addItem("VtxChisq",     &VTX_Chisq);

	addItem("VtxS1",        &VTX_S1);       
	addItem("VtxS2",        &VTX_S2);       
	addItem("VtxAddedRL",   &VTX_AddedRL); 

    addItem("Vtx2TransDoca", &VTX2_transDoca);
    addItem("Vtx2LongDoca",  &VTX2_longDoca);

	addItem("Vtx2XDir",     &VTX2_xdir);     
	addItem("Vtx2YDir",     &VTX2_ydir);     
	addItem("Vtx2ZDir",     &VTX2_zdir);     
    addItem("Vtx2X0",       &VTX2_x0);       
	addItem("Vtx2Y0",       &VTX2_y0);       
	addItem("Vtx2Z0",       &VTX2_z0);

    addItem("Vtx2Angle",     &VTX2_Angle);    
	addItem("Vtx2DOCA",      &VTX2_DOCA);
	addItem("Vtx2HeadSep",   &VTX2_Head_Sep);
    addItem("Vtx2Status",    &VTX2_Status);

	addItem("VtxNeutXDir" , &VTXN_xdir);
	addItem("VtxNeutYDir" , &VTXN_ydir);
	addItem("VtxNeutZDir" , &VTXN_zdir);
	addItem("VtxNeutSXX",   &VTXN_Sxx);
    addItem("VtxNeutSXY",   &VTXN_Sxy);
    addItem("VtxNeutSYY",   &VTXN_Syy);
    addItem("VtxNeutChgWt",   &VTXN_ChgWt);
    addItem("VtxNeutNeutWt",   &VTXN_NeutWt);

	addItem("VtxNeut1XDir" , &VTXN1_xdir);
	addItem("VtxNeut1YDir" , &VTXN1_ydir);
	addItem("VtxNeut1ZDir" , &VTXN1_zdir);
	addItem("VtxNeut1SXX",   &VTXN1_Sxx);
    addItem("VtxNeut1SXY",   &VTXN1_Sxy);
    addItem("VtxNeut1SYY",   &VTXN1_Syy);
    addItem("VtxNeut1ChgWt",   &VTXN1_ChgWt);
    addItem("VtxNeut1NeutWt",   &VTXN1_NeutWt);

	zeroVals();

	return sc;
}


StatusCode VtxValsTool::calculate()
{
	StatusCode sc = StatusCode::SUCCESS;

	// Recover Track associated info. 
	SmartDataPtr<Event::TkrTrackCol>  
		pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol);
	SmartDataPtr<Event::TkrVertexCol>     
		pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);

	if(!pVerts) return sc; 

    VTX_numVertices = (int) pVerts->size();

	// Get the first Vertex - First track of first vertex = Best Track
	Event::TkrVertexConPtr pVtxr = pVerts->begin(); 
	if(pVtxr == pVerts->end()) return sc; 

	Event::TkrVertex*   gamma = *pVtxr++; 
	SmartRefVector<Event::TkrTrack>::const_iterator pTrack1 = gamma->getTrackIterBegin(); 
	const Event::TkrTrack* track_1 = *pTrack1;

	int nParticles = gamma->getNumTracks(); 

	Point  x0 = gamma->getPosition();
	Vector t0 = gamma->getDirection();

	VTX_xdir      = t0.x();
	VTX_ydir      = t0.y();
	VTX_zdir      = t0.z();

	VTX_Phi       = (-t0).phi();
	if (VTX_Phi<0.0f) VTX_Phi += static_cast<float>(2*M_PI);
	VTX_Theta     = (-t0).theta();

    const Event::TkrTrackParams& VTX_Cov = gamma->getVertexParams();
    VTX_Sxx         = VTX_Cov.getxSlpxSlp();
    VTX_Sxy         = VTX_Cov.getxSlpySlp();
    VTX_Syy         = VTX_Cov.getySlpySlp();
    double sinPhi     = sin(VTX_Phi);
    double cosPhi     = cos(VTX_Phi);
    VTX_ThetaErr      = t0.z()*t0.z()*sqrt(std::max(0.0, cosPhi*cosPhi*VTX_Sxx + 
        2.*sinPhi*cosPhi*VTX_Sxy + sinPhi*sinPhi*VTX_Syy)); 
    VTX_PhiErr        = (-t0.z())*sqrt(std::max(0.0, sinPhi*sinPhi*VTX_Sxx - 
        2.*sinPhi*cosPhi*VTX_Sxy + cosPhi*cosPhi*VTX_Syy));
    //VTX_ErrAsym     = fabs(VTX_Sxy/(VTX_Sxx + VTX_Syy));
    //VTX_CovDet      = sqrt(std::max(0.0f,VTX_Sxx*VTX_Syy-VTX_Sxy*VTX_Sxy))*VTX_zdir*VTX_zdir;



	VTX_x0        = x0.x();
	VTX_y0        = x0.y();
	VTX_z0        = x0.z();

    VTX_Status  = gamma->getStatusBits(); 

	// Get the first track location and direction
	Point  x1 = track_1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
	Vector t1 = track_1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);

	// Check if there is a second track in the event.  This track may not be 
	// associated with the first track to form the first vertex.

    double cost1t2, t1t2;
    Point x2H;

	if(nParticles > 1) {
		pTrack1++;
        const Event::TkrTrack* track_2 = *pTrack1;

		Point  x2 = track_2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
	    Vector t2 = track_2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);

		x2H = x2 + ((x1.z()-x2.z())/t2.z()) * t2;

		VTX_Head_Sep = (x1-x2H).magnitude(); 

		double cost1t2 = t1*t2; 
		double t1t2  = acos(cost1t2); 
		VTX_Angle = t1t2;
		VTX_DOCA  = gamma->getDOCA(); 
		VTX_S1    = gamma->getTkr1ArcLen();
		VTX_S2    = gamma->getTkr2ArcLen();

		// Set a rogue value here in case this is a single 
		if(VTX_xdir == t1.x() && VTX_ydir == t1.y()) VTX_Angle = -.1f;

		VTX_Quality = gamma->getQuality(); 
		VTX_Chisq   = gamma->getChiSquare(); 
		VTX_AddedRL = gamma->getAddedRadLen();
	}

    if(pVerts->size()>1) {
        Event::TkrVertex*   vtx2 = *pVtxr++; 

		if(!(vtx2->getStatusBits()& Event::TkrVertex::NEUTRALVTX)) {

	    Point  x2 = vtx2->getPosition();
	    Vector t2 = vtx2->getDirection();
        VTX2_Status    = vtx2->getStatusBits();

	    VTX2_xdir      = t2.x();
	    VTX2_ydir      = t2.y();
	    VTX2_zdir      = t2.z();

	    VTX2_x0        = x2.x();
	    VTX2_y0        = x2.y();
	    VTX2_z0        = x2.z();

        Doca vtx0(x0, t0);
        VTX2_transDoca = vtx0.docaOfPoint(x2);
        VTX2_longDoca  = vtx0.arcLenRay1();

        pTrack1 = vtx2->getTrackIterBegin(); 
        nParticles = vtx2->getNumTracks(); 
        const Event::TkrTrack* track_1a = *pTrack1;
        // Get the first track location and direction
        Point  x1 = track_1a->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
        Vector t1 = track_1a->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
        if(nParticles > 1) {
            pTrack1++;
            const Event::TkrTrack* track_2a = *pTrack1;

            x2 = track_2a->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
            t2 = track_2a->front()->getDirection(Event::TkrTrackHit::SMOOTHED);

            x2H = x2 + ((x1.z()-x2.z())/t2.z()) * t2;

            VTX2_Head_Sep = (x1-x2H).magnitude(); 

            cost1t2 = t1*t2; 
            t1t2  = acos(cost1t2); 
            VTX2_Angle = t1t2;
            VTX2_DOCA  = vtx2->getDOCA();

 		    // Set a rogue value here in case this is a single 
		    if(VTX2_xdir == t1.x() && VTX2_ydir == t1.y()) VTX2_Angle = -.1f;
       }
		}}

	//Neutral Vertex section
   	Event::TkrVertexConPtr pVtxN = pVerts->begin(); 
	Event::TkrVertex*   vtxN; 
	bool VTX_Set = false; 
	for(; pVtxN != pVerts->end(); pVtxN++)
    {	    
        vtxN = *pVtxN; 
		const Event::TkrTrackParams& VTXN_Cov = vtxN->getVertexParams();
		if(vtxN->getStatusBits()& Event::TkrVertex::NEUTRALVTX) {
			Vector tN = vtxN->getDirection();
			float chrg_wt = vtxN->getTkr1ArcLen();
			float neut_wt = vtxN->getTkr2ArcLen();

			if(vtxN->getStatusBits()& Event::TkrVertex::ONETKRVTX){
				if(!VTX_Set) {
					VTXN_xdir      = tN.x();
					VTXN_ydir      = tN.y();
					VTXN_zdir      = tN.z();
					VTXN_Sxx       = VTXN_Cov.getxSlpxSlp();
					VTXN_Sxy       = VTXN_Cov.getxSlpySlp();
					VTXN_Syy       = VTXN_Cov.getySlpySlp();
					VTXN_ChgWt     = chrg_wt;
					VTXN_NeutWt    = neut_wt;
					VTX_Set = true;
				}
				VTXN1_xdir      = tN.x();
				VTXN1_ydir      = tN.y();
				VTXN1_zdir      = tN.z();
				VTXN1_Sxx       = VTXN_Cov.getxSlpxSlp();
				VTXN1_Sxy       = VTXN_Cov.getxSlpySlp();
				VTXN1_Syy       = VTXN_Cov.getySlpySlp();
				VTXN1_ChgWt     = chrg_wt;
				VTXN1_NeutWt    = neut_wt;
			} else {
				VTXN_xdir      = tN.x();
				VTXN_ydir      = tN.y();
				VTXN_zdir      = tN.z();
				VTXN_Sxx       = VTXN_Cov.getxSlpxSlp();
				VTXN_Sxy       = VTXN_Cov.getxSlpySlp();
				VTXN_Syy       = VTXN_Cov.getySlpySlp();
				VTXN_ChgWt     = chrg_wt;
				VTXN_NeutWt    = neut_wt;
				VTX_Set = true;
	}   }	}

	return sc;
}
