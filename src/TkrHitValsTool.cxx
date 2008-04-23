/** @file TkrHitValsTool.cxx
@brief Calculates the Tkr hit analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrHitValsTool.cxx,v 1.14 2008/04/21 21:05:35 lsrea Exp $
*/

// Include files

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "TkrUtil/ITkrQueryClustersTool.h"

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
    int Tkr_Cnv_Lyr_Hits;
    int Tkr_Max_controller_hits;
    int Tkr_Fst_Cnv_Lyr;
    int Tkr_NCnv_Lyrs_Hit;
    
    int Tkr_HitsPerLyr[_nLayers];

    ITkrQueryClustersTool* m_clusTool;
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

/** @page anatup_vars 
    @section tkrhitvalstool TkrHitValsTool Variables

<table>
<tr><th> Variable <th> Type  <th> Description					
<tr><td> TkrNumHits 	
<td>I<td>   Total number of TKR clusters 
<tr><td> TkrFirstLayer
<td>I<td>   First layer containing a cluster 
<tr><td> TkrNumLayersHit
<td>I<td>   Total number of hit layers 
<tr><td> TkrHitsInLyrNN, NN=(00,17)   
<td>I<td>   Number of clusters in (bi)layer NN 
           (numbered from the bottom of the tracker) 
</table>
*/

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

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }
    
    // load up the map

    addItem("TkrNumHits",            &Tkr_Cnv_Lyr_Hits);       
    addItem("TkrFirstLayer",         &Tkr_Fst_Cnv_Lyr);        
    addItem("TkrNumLayersHit",       &Tkr_NCnv_Lyrs_Hit);

    int i;
    char buffer[20];
    for(i=0;i<_nLayers;++i) {
        // vars of form TkrHitsInLyrNN, NN = (00,17), for _nLayers==18
        sprintf(buffer, "TkrHitsInLyr%02i",i);
        addItem(buffer, &Tkr_HitsPerLyr[i]);
    }
            
    zeroVals();
    
    return sc;
}


StatusCode TkrHitValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    if (!pClusters) return sc;

    //Make sure we have valid cluster data
    if (pClusters)
    {
        int layerIdx = _nLayers;
        while(layerIdx--)
        {
            int hitCount = m_clusTool->getClusters(idents::TkrId::eMeasureX,layerIdx).size()
                         + m_clusTool->getClusters(idents::TkrId::eMeasureY,layerIdx).size();
            
            if (hitCount > 0)
            {
                Tkr_Fst_Cnv_Lyr    = layerIdx;
                Tkr_NCnv_Lyrs_Hit += 1;
            }
            
            Tkr_HitsPerLyr[layerIdx] = hitCount;
        }
        
        Tkr_Cnv_Lyr_Hits = pClusters->size();
    }
    
    return sc;
}
