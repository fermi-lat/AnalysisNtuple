// $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrHitValsTool.cxx,v 1.2 2003/02/24 17:55:45 lsrea Exp $

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

    m_ntupleMap["TKR_Cnv_Lyr_Hits"]        = &Tkr_Cnv_Lyr_Hits;       
    m_ntupleMap["TKR_Max_controller_hits"] = &Tkr_Max_controller_hits;
    m_ntupleMap["TKR_Fst_Cnv_Lyr"]         = &Tkr_Fst_Cnv_Lyr;        
    m_ntupleMap["TKR_NCnv_Lyrs_Hit"]       = &Tkr_NCnv_Lyrs_Hit;      
    
    m_ntupleMap["TKR_Hits_In_Lyr_0"] =      &Tkr_HitsPerLyr[0];      
    m_ntupleMap["TKR_Hits_In_Lyr_1"] =      &Tkr_HitsPerLyr[1];      
    m_ntupleMap["TKR_Hits_In_Lyr_2"] =      &Tkr_HitsPerLyr[2];    
    m_ntupleMap["TKR_Hits_In_Lyr_3"] =      &Tkr_HitsPerLyr[3];      
    m_ntupleMap["TKR_Hits_In_Lyr_4"] =      &Tkr_HitsPerLyr[4];      
    m_ntupleMap["TKR_Hits_In_Lyr_5"] =      &Tkr_HitsPerLyr[5];      
    m_ntupleMap["TKR_Hits_In_Lyr_6"] =      &Tkr_HitsPerLyr[6];      
    m_ntupleMap["TKR_Hits_In_Lyr_7"] =      &Tkr_HitsPerLyr[7];      
    m_ntupleMap["TKR_Hits_In_Lyr_8"] =      &Tkr_HitsPerLyr[8];      
    m_ntupleMap["TKR_Hits_In_Lyr_9"] =      &Tkr_HitsPerLyr[9];      
    m_ntupleMap["TKR_Hits_In_Lyr_10"] =     &Tkr_HitsPerLyr[10];     
    m_ntupleMap["TKR_Hits_In_Lyr_11"] =     &Tkr_HitsPerLyr[11];     
    m_ntupleMap["TKR_Hits_In_Lyr_12"] =     &Tkr_HitsPerLyr[12];     
    m_ntupleMap["TKR_Hits_In_Lyr_13"] =     &Tkr_HitsPerLyr[13];     
    m_ntupleMap["TKR_Hits_In_Lyr_14"] =     &Tkr_HitsPerLyr[14];     
    m_ntupleMap["TKR_Hits_In_Lyr_15"] =     &Tkr_HitsPerLyr[15];     
    m_ntupleMap["TKR_Hits_In_Lyr_16"] =     &Tkr_HitsPerLyr[16];     
    m_ntupleMap["TKR_Hits_In_Lyr_17"] =     &Tkr_HitsPerLyr[17]; 
    
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
