
/** @file VtxValsTool.cxx
@brief Calculates the Vtx analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/VtxValsTool.cxx,v 1.14 2004/09/30 17:12:21 lsrea Exp $
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
	double VTX_xdir;
	double VTX_ydir;
	double VTX_zdir;
	double VTX_Phi;
	double VTX_Theta;
	double VTX_x0;
	double VTX_y0;
	double VTX_z0;
	double VTX_Angle;
	double VTX_DOCA;
	double VTX_Head_Sep;

	double VTX_S1;
	double VTX_S2;

	double VTX_Quality; 
	double VTX_Chisq; 
	double VTX_AddedRL;
	double VTX_Status;

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
	if (VTX_Phi<0.0) VTX_Phi += 2*M_PI;
	VTX_Theta     = (-t0).theta();

	VTX_x0        = x0.x();
	VTX_y0        = x0.y();
	VTX_z0        = x0.z();

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

		VTX_Head_Sep = (x1-x2).magnitude(); 

		double cost1t2 = t1*t2; 
		double t1t2  = acos(cost1t2); 
		VTX_Angle = t1t2;
		VTX_DOCA  = gamma->getDOCA(); 
		VTX_S1    = gamma->getTkr1ArcLen();
		VTX_S2    = gamma->getTkr2ArcLen();

		// Set a rogue value here in case this is a single 
		if(VTX_xdir == t1.x() && VTX_ydir == t1.y()) VTX_Angle = -.1;

		VTX_Quality = gamma->getQuality(); 
		VTX_Chisq   = gamma->getChiSquare(); 
		VTX_Status  = gamma->getStatusBits(); 
		VTX_AddedRL = gamma->getAddedRadLen();
	}

	return sc;
}
