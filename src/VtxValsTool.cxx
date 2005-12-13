
/** @file VtxValsTool.cxx
@brief Calculates the Vtx analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/VtxValsTool.cxx,v 1.20 2005/11/18 16:49:15 lsrea Exp $
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
	float VTX_xdir;
	float VTX_ydir;
	float VTX_zdir;
	float VTX_Phi;
	float VTX_Theta;
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

};

// Static factory for instantiation of algtool objects
static ToolFactory<VtxValsTool> s_factory;
const IToolFactory& VtxValsToolFactory = s_factory;

// Standard Constructor
VtxValsTool::VtxValsTool(const std::string& type, 
						 const std::string& name, 
						 const IInterface* parent)
						 : ValBase( type, name, parent )
{    
	// Declare additional interface
	declareInterface<IValsTool>(this); 
}

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

    /** @page anatup_vars 
    @section vtxvalstool VtxValsTool Variables

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> Vtx[X/Y/Z]Dir 
<td> F
<td>        [x/y/z] direction cosine of the (1st) vertex 
<tr><td> VtxPhi 
<td> F
<td>        Azimuthal angle of vertex, radians 
            (direction of source, not flight direction!) Range: (0,2pi) 
<tr><td> VtxTheta 
<td> F
<td>        Polar angle of vertex, radians (ditto direction)
<tr><td> Vtx[X/Y/Z]0 
<td> F
<td>        [x/y/z] coordinate of vertex; if the two tracks making up 
            the vertex are nearly parallel, 
            the coordinates of the vertex may become very large. 
<tr><td> VtxAngle 
<td> F
<td>        Angle between the two tracks of the vertex (radians) 
<tr><td> VtxDOCA 
<td> F
<td>        Distance of closest approach between the two tracks 
<tr><td> VtxHeadSep 
<td> F
<td>        Distance between the heads of the two tracks 
<tr><td> VtxStatus 
<td> F
<td>        Summary of track composition and topology.
            See TkrVertex.h in the Event package for the current description.
            The definitions as of GR v7r2 are:
@verbatim

    |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
                                         [Track Topology ] [Track composition]

    ONETKRVTX = 0x0001,  //Set if single track vertex
    TWOTKRVTX = 0x0002,  //Set if 2 track vertex
    MULTKRVTX = 0x0004,  //Set if >2 track vertex
    DOCAVTX   = 0x0010,  //Set if vertex location set by DOCA point
    FIRSTHIT  = 0x0020,  //Set if two tracks share first hit
    STAGVTX   = 0x0040,  //Set if tracks don't start in same plane (staggered)
    CROSSTKR  = 0x0080   //Set if DOCA location lies inside track hits

@endverbatim
<tr><td> VtxQuality 
<td> F
<td>        Vertex quality parameter used to order the possible vertices and 
            select the best one. <strong>Should generally not be used in analysis</strong>. 
<tr><td> VtxChisq 
<td> F
<td>         The covariant chi-squared for the pairing of the tracks. 
<tr><td> VtxS[1/2] 
<td> F
<td>        Distance of DOCA point from head of track [1/2] 
<tr><td> VtxAddedRL 
<td> F
<td>        The additional radiation lengths prior to the first measured silicon 
            strip hit at the vertex location. <em>New!</em> 
</table>
    */

	// Pair reconstruction
	addItem("VtxXDir",      &VTX_xdir);     
	addItem("VtxYDir",      &VTX_ydir);     
	addItem("VtxZDir",      &VTX_zdir);     
	addItem("VtxPhi",       &VTX_Phi);  
	addItem("VtxTheta",     &VTX_Theta);  
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

	VTX_x0        = x0.x();
	VTX_y0        = x0.y();
	VTX_z0        = x0.z();

    VTX_Status  = gamma->getStatusBits(); 

	// Get the first track location and direction
	Point  x1 = track_1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
	Vector t1 = track_1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);

	// Check if there is a second track in the event.  This track may not be 
	// associated with the first track to form the first vertex.
	if(nParticles > 1) {
		pTrack1++;
		const Event::TkrTrack* track_2 = *pTrack1;

		Point  x2 = track_2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
	    Vector t2 = track_2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);

		Point x2H = x2 + ((x1.z()-x2.z())/t2.z()) * t2;

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
    
	return sc;
}
