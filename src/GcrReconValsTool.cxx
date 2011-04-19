    
/** @file GcrReconValsTool.cxx
@brief Calculates the GcrRecon analysis variables
@author C. Lavalley

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/GcrReconValsTool.cxx,v 1.1.618.1 2010/09/18 02:35:14 heather Exp $
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

//@@@CL 06/26/06
#include "Event/Recon/CalRecon/GcrReconClasses.h"
//@@@CL 06/26/06

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

/*! @class GcrReconValsTool
@brief calculates Cal Mip values

  @author C.Lavalley
  */
  
class GcrReconValsTool :   public ValBase
{
public:  
    GcrReconValsTool( const std::string& type, 
                    const std::string& name, 
                    const IInterface* parent);
      
    virtual ~GcrReconValsTool() { }
      
    StatusCode initialize();
      
    StatusCode calculate();
      
private:

    static const int NTOW = 16;
    static const int NLAY = 8;
    static const int NCOL = 12;

    double m_gcrRecon[NTOW*NLAY*NCOL]; 
    float m_inferedZ;
    
};
  
// Static factory for instantiation of algtool objects
//static ToolFactory<GcrReconValsTool> s_factory;
//const IToolFactory& GcrReconValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(GcrReconValsTool);
  
// Standard Constructor
GcrReconValsTool::GcrReconValsTool(const std::string& type, 
                               const std::string& name, 
                               const IInterface* parent)
               : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}
  
StatusCode GcrReconValsTool::initialize()
{


//std::cout << "BEGIN initialize in GcrReconValsTool" << std::endl;
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;


    
    addItem("GcrRecon[1536]",     m_gcrRecon);
//std::cout << "END initialize in GcrReconValsTool" << std::endl;
 
    return sc;
}

StatusCode GcrReconValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the calMipTrack collection   

    MsgStream log(msgSvc(), name());
    //log << MSG::INFO << "BEGIN calculate in GcrReconValsTool" << endreq;
    
    //INITIALISATION OF m_gcrRecon
    int j;
   for(int itow=0; itow<NTOW;itow++) 
      for (int ilay=0;ilay<NLAY;ilay++)
        for(int icol=0;icol<NCOL;icol++)
	{
    	  j=(itow*NLAY+ilay)*NCOL+icol;
	  m_gcrRecon[j]=-1000;
    
        }
	
    
    SmartDataPtr<Event::GcrXtalCol> p_gcrXtalCol(m_pEventSvc, EventModel::CalRecon::GcrXtalCol); 
    m_inferedZ=-1000;
    if(p_gcrXtalCol){
	    int i=0;
	    int itow,ilay,icol;
	    for(Event::GcrXtalCol::const_iterator gcrXtalIter=p_gcrXtalCol->begin(); gcrXtalIter != p_gcrXtalCol->end(); gcrXtalIter++)
	    {
	      Event::GcrXtal* p_gcrXtal = *gcrXtalIter;
	      
	      //Event::CalXtalRecData* xtalData = p_gcrSelectedXtal->getXtal();
	      idents::CalXtalId xtalId = p_gcrXtal->getXtalId();
	      itow = xtalId.getTower();
	      ilay = xtalId.getLayer();
	      icol = xtalId.getColumn();

	      i=(itow*NLAY+ilay)*NCOL+icol;
	      //log << MSG::INFO << "GcrSelValsTool::calculate p_gcrSelectedXtal->getInferedZ()=" <<  p_gcrSelectedXtal->getInferedZ()<< endreq;
	      
	      /**std::printf("%10.4f \n", p_gcrSelectedXtal->getPathLength());
	      fflush(stdout);*/
	      
	      m_gcrRecon[i] = p_gcrXtal->getPathLength();
          
	    }

    }
    else{
      log << MSG::INFO << "no gcrXtalCol found " << endreq;
      
      }
    /**FOR DEBUG: DISPLAY m_gcrRecon
    int i;
    for(int itow=0; itow<NTOW;itow++) 
      for (int ilay=0;ilay<NLAY;ilay++)
        for(int icol=0;icol<NCOL;icol++)
	{
    	  i=(itow*NLAY+ilay)*NCOL+icol;
	  if(m_gcrRecon[i]>0)
	    log << MSG::INFO << "m_gcrRecon["<<itow<<","<<icol<<","<<ilay <<"]= "<< m_gcrRecon[i] << endreq;
    
        }*/
      
    //log << MSG::INFO << "END calculate in GcrReconValsTool End  " << endreq;
    return sc;
}


/**

Mail from Leon


>
> AnaTup is already set up to write out arrays... the part that needs some work is to access the array elements from *inside* Gleam.
>
> But I assume you don't need to do that!
>
> So here's how it works:
>
>
> float myVal[100];    // or double, int
>
> ...
>
> AddItem("MyVal[100]", myVal); // note that myVal is a pointer, could also be &myVal[0]
>
> ...
>
> int i;
> for(i=0; i<100; ++i) {
>     myVal[i] = i*3.14159;
> }
>
> ...
>
> Then MyVal[100] should be available in the ntuple.
>
> L.
>
>
> 

*/
