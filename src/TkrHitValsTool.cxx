// $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrHitValsTool.cxx,v 1.2 2003/02/26 17:38:48 lsrea Exp $

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

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"

namespace {
    const int _nLayers = 18;
}

/*! @class TkrHitValsTool
@brief calculates TkrHit values

@authors Bill Atwood, Leon Rochester
*/

class TkrHitValsTool : public ValBase
{
public:
    
    TkrHitValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~TkrHitValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    
    //TkrClusters Tuple Items
    double Tkr_Cnv_Lyr_Hits;
    double Tkr_Max_controller_hits;
    double Tkr_Fst_Cnv_Lyr;
    double Tkr_NCnv_Lyrs_Hit;
    
    double Tkr_HitsPerLyr[_nLayers];
};

// Static factory for instantiation of algtool objects
static ToolFactory<TkrHitValsTool> s_factory;
const IToolFactory& TkrHitValsToolFactory = s_factory;

// Standard Constructor
TkrHitValsTool::TkrHitValsTool(const std::string& type, 
                               const std::string& name, 
                               const IInterface* parent)
                               : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode TkrHitValsTool::initialize()
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

    addItem("TKR_Cnv_Lyr_Hits",        &Tkr_Cnv_Lyr_Hits);       
    addItem("TKR_Max_controller_hits", &Tkr_Max_controller_hits);
    addItem("TKR_Fst_Cnv_Lyr",         &Tkr_Fst_Cnv_Lyr);        
    addItem("TKR_NCnv_Lyrs_Hit",       &Tkr_NCnv_Lyrs_Hit);      
    
    addItem("TKR_Hits_In_Lyr_0",      &Tkr_HitsPerLyr[0]);      
    addItem("TKR_Hits_In_Lyr_1",      &Tkr_HitsPerLyr[1]);      
    addItem("TKR_Hits_In_Lyr_2",      &Tkr_HitsPerLyr[2]);    
    addItem("TKR_Hits_In_Lyr_3",      &Tkr_HitsPerLyr[3]);      
    addItem("TKR_Hits_In_Lyr_4",      &Tkr_HitsPerLyr[4]);      
    addItem("TKR_Hits_In_Lyr_5",      &Tkr_HitsPerLyr[5]);      
    addItem("TKR_Hits_In_Lyr_6",      &Tkr_HitsPerLyr[6]);      
    addItem("TKR_Hits_In_Lyr_7",      &Tkr_HitsPerLyr[7]);      
    addItem("TKR_Hits_In_Lyr_8",      &Tkr_HitsPerLyr[8]);      
    addItem("TKR_Hits_In_Lyr_9",      &Tkr_HitsPerLyr[9]);      
    addItem("TKR_Hits_In_Lyr_10",     &Tkr_HitsPerLyr[10]);     
    addItem("TKR_Hits_In_Lyr_11",     &Tkr_HitsPerLyr[11]);     
    addItem("TKR_Hits_In_Lyr_12",     &Tkr_HitsPerLyr[12]);     
    addItem("TKR_Hits_In_Lyr_13",     &Tkr_HitsPerLyr[13]);     
    addItem("TKR_Hits_In_Lyr_14",     &Tkr_HitsPerLyr[14]);     
    addItem("TKR_Hits_In_Lyr_15",     &Tkr_HitsPerLyr[15]);     
    addItem("TKR_Hits_In_Lyr_16",     &Tkr_HitsPerLyr[16]);     
    addItem("TKR_Hits_In_Lyr_17",     &Tkr_HitsPerLyr[17]); 
    
    zeroVals();
    
    return sc;
}


StatusCode TkrHitValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>    
        pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    
    //Make sure we have valid cluster data
    if (pClusters)
    {
        int layerIdx = _nLayers;
        while(layerIdx--)
        {
            int hitCount = pClusters->nHits(Event::TkrCluster::X, layerIdx) 
                + pClusters->nHits(Event::TkrCluster::Y, layerIdx);
            
            if (hitCount > 0)
            {
                Tkr_Fst_Cnv_Lyr    = layerIdx;
                Tkr_NCnv_Lyrs_Hit += 1;
            }
            
            Tkr_HitsPerLyr[layerIdx] = hitCount;
        }
        
        Tkr_Cnv_Lyr_Hits = pClusters->nHits();
    }
    
    return sc;
}
