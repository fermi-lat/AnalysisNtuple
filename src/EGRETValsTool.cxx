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

#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

// Allow access to EGRET vars
#include "Event/EGRET/EGRETGamma.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

class EGRETValsTool : public ValBase
{
public:

	EGRETValsTool( const std::string& type, 
		const std::string& name, 
		const IInterface* parent);

	virtual ~EGRETValsTool() { }

	StatusCode initialize();

	StatusCode calculate();

private:

	//Pure EGRET Tuple Items
	double EGRET_Id;
	double EGRET_Energy;
	double EGRET_EFrac;
	double EGRET_OpenAngle; 


	double EGRET_x0;
	double EGRET_y0;
	double EGRET_z0;

	double EGRET_xdir;
	double EGRET_ydir;
	double EGRET_zdir;

	//EGRET - Compared to Recon Items
	double EGRET_x_err; 
	double EGRET_y_err;
	double EGRET_z_err;

	double EGRET_xdir_err; 
	double EGRET_ydir_err;
	double EGRET_zdir_err;

	double EGRET_dir_err;
	double EGRET_TKR1_dir_err;
	double EGRET_TKR2_dir_err;

//	char*  EGRET_Pulsar;
	double EGRET_PulsarPhase;


	// to decode the particle charge
	//  IParticlePropertySvc* m_ppsvc;    
};

// Static factory for instantiation of algtool objects
static ToolFactory<EGRETValsTool> s_factory;
const IToolFactory& EGRETValsToolFactory = s_factory;

// Standard Constructor
EGRETValsTool::EGRETValsTool(const std::string& type, 
							 const std::string& name, 
							 const IInterface* parent)
							 : ValBase( type, name, parent )
{    
	// Declare additional interface
	declareInterface<IValsTool>(this); 
}

StatusCode EGRETValsTool::initialize()
{
	StatusCode sc = StatusCode::SUCCESS;

	MsgStream log(msgSvc(), name());

	if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

	// get the services    
	/*
	if( serviceLocator() ) {
	if( service("ParticlePropertySvc", m_ppsvc, true).isFailure() ) {
	log << MSG::ERROR << "Service [ParticlePropertySvc] not found" << endreq;
	}
	} else {
	return StatusCode::FAILURE;
	}
	*/
	// load up the map

	addItem("EGRETEnergy",       &EGRET_Energy);  
	addItem("EGRETEFrac",        &EGRET_EFrac);
	addItem("EGRETOpenAngle",    &EGRET_OpenAngle);
	addItem("EGRETX0",           &EGRET_x0);           
	addItem("EGRETY0",           &EGRET_y0);           
	addItem("EGRETZ0",           &EGRET_z0);           

	addItem("EGRETXDir",         &EGRET_xdir);         
	addItem("EGRETYDir",         &EGRET_ydir);         
	addItem("EGRETZDir",         &EGRET_zdir);         

	addItem("EGRETXErr",         &EGRET_x_err);        
	addItem("EGRETYErr",         &EGRET_y_err);        
	addItem("EGRETZErr",         &EGRET_z_err);        

	addItem("EGRETXDirErr",      &EGRET_xdir_err);     
	addItem("EGRETYDirErr",      &EGRET_ydir_err);     
	addItem("EGRETZDirErr",      &EGRET_zdir_err);     

	addItem("EGRETDirErr",       &EGRET_dir_err);      
	addItem("EGRETTkr1DirErr",   &EGRET_TKR1_dir_err); 
	addItem("EGRETTkr2DirErr",   &EGRET_TKR2_dir_err);

//	addItem("EGRETPulsar",       &EGRET_Pulsar);
	addItem("EGRETPulsarPhase",  &EGRET_PulsarPhase);

	zeroVals();

	return sc;
}


StatusCode EGRETValsTool::calculate()
{
	StatusCode sc = StatusCode::SUCCESS;

	// Recover EGRET Track associated info. 
	SmartDataPtr<Event::TkrPatCandCol>    eTracks(m_pEventSvc,EventModel::TkrRecon::TkrPatCandCol);
	SmartDataPtr<Event::EGRETGammaCol>    eVerts(m_pEventSvc,EventModel::Digi::EGRETGammaCol);

	// Recover Gleam Track associated info. 
	SmartDataPtr<Event::TkrFitTrackCol> 
		pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
	SmartDataPtr<Event::TkrVertexCol>     
		pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);

	if (eVerts) 
	{
		// Position of Gamma start
		Event::EGRETGamma* egamma = *(eVerts->begin());
		Point Egret_x0 = egamma->getPosition();

		Vector Egret_t0 = egamma->getDirection();



		//Pure EGRET Tuple Items
		EGRET_Energy = egamma->getEnergy();

		EGRET_PulsarPhase = egamma->getPulsarPhase();

		EGRET_x0     = Egret_x0.x();
		EGRET_y0     = Egret_x0.y();
		EGRET_z0     = Egret_x0.z();

		EGRET_xdir   = Egret_t0.x();
		EGRET_ydir   = Egret_t0.y();
		EGRET_zdir   = Egret_t0.z();

		EGRET_OpenAngle = egamma->getOpeningAngle();

		if(eTracks) {

			Event::TkrPatCand* etkr1 = *(eTracks->begin());
			double e1 = etkr1->getEnergy();
			EGRET_EFrac = e1/EGRET_Energy;
			if(eTracks->size() >1) 
			{
				Event::TkrPatCand* etkr2 = *(eTracks->begin()+1);		
				double e2 = etkr2->getEnergy();

				if(e1 < e2) EGRET_EFrac = e2/EGRET_Energy;
			}

		}

		//Make sure we have valid reconstructed data
		if (pVerts->size()>0) {

			// Get the first Vertex - First track of first vertex = Best Track
			Event::TkrVertexConPtr pVtxr = pVerts->begin();
			Event::TkrVertex*   gamma = *(pVtxr);
			Point  x0 = gamma->getPosition();
			Vector t0 = gamma->getDirection();


			// Reference position errors at the start of recon track(s)
			double arc_len = (x0.z()-Egret_x0.z())/Egret_t0.z();
			Point x_start = Egret_x0 + arc_len*Egret_t0;
			EGRET_x_err  = x0.x()-x_start.x(); 
			EGRET_y_err  = x0.y()-x_start.y();
			// except for z, use the difference between EGRET conversion point and start of vertex track
			EGRET_z_err  = x0.z()-Egret_x0.z();

			EGRET_xdir_err = t0.x()-Egret_t0.x(); 
			EGRET_ydir_err = t0.y()-Egret_t0.y();
			EGRET_zdir_err = t0.z()-Egret_t0.z();

			double cost0tEGRET = t0*Egret_t0;

			EGRET_dir_err  = acos(cost0tEGRET);

			int nParticles = gamma->getNumTracks(); 

			SmartRefVector<Event::TkrFitTrackBase>::const_iterator pTrack1 = gamma->getTrackIterBegin();  
			const Event::TkrFitTrackBase* trackBase = *pTrack1;
			SmartRef<Event::TkrKalFitTrack> track_1  = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase); 

			Point  x1 = track_1->getPosition();
			Vector t1 = track_1->getDirection();
			double cost1tEGRET = t1*Egret_t0;

			EGRET_TKR1_dir_err  = acos(cost1tEGRET);

			if(nParticles > 1) {
				pTrack1++;
				trackBase = *pTrack1;
				SmartRef<Event::TkrKalFitTrack> track_2 = dynamic_cast<const Event::TkrKalFitTrack*>(trackBase);
				Point  x2 = track_2->getPosition();
				Vector t2 = track_2->getDirection();
				double cost2tEGRET = t2*Egret_t0;

				EGRET_TKR2_dir_err  = acos(cost2tEGRET);
			}

		} 
	}

	return sc;
}
