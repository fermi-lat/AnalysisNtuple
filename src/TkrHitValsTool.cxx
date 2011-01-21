/** @file TkrHitValsTool.cxx
@brief Calculates the Tkr hit analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TkrHitValsTool.cxx,v 1.22 2010/12/06 20:59:14 lsrea Exp $
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
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"

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
    int Tkr_numHitsOnTracks;
    int Tkr_numGhosts;
    int Tkr_numDiags;
    int Tkr_numBoth;
    int Tkr_numToT255s;
    int Tkr_numSaturated;
    int Tkr_numWideClusters;
    int Tkr_numWiderClusters;
    int Tkr_numSaturatedGhosts;
    int Tkr_numWideGhostClusters;
    int Tkr_numWiderGhostClusters;
    int Tkr_numGhostsOnTracks;
    int Tkr_numDiagsOnTracks;
    int Tkr_numBothOnTracks;
    int Tkr_numToT255sOnTracks;
    int Tkr_numFlaggedTrackHits;
    int Tkr_numSaturatedOnTracks;
    int Tkr_numWideClustersOnTracks;
    int Tkr_numWiderClustersOnTracks;
    int Tkr_numSaturatedGhostsOnTracks;
    int Tkr_numWideGhostClustersOnTracks;
    int Tkr_numWiderGhostClustersOnTracks;
    int Tkr_Max_controller_hits;
    int Tkr_Fst_Cnv_Lyr;
    int Tkr_NCnv_Lyrs_Hit;

    int Tkr_HitsPerLyr[_nLayers];
    int Tkr_numRCTruncs;
    int Tkr_numCCTruncs;

    ITkrQueryClustersTool* m_clusTool;

    int m_minWide;
    int m_minWider;
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

    declareProperty("minWide", m_minWide=5);
    declareProperty("minWider", m_minWider=9);

}

/** @page anatup_vars 
@section tkrhitvalstool TkrHitValsTool Variables

<table>
<tr><th> Variable <th> Type  <th> Description                                        
<tr><td> TkrNumHits         
<td>I<td>   Total number of TKR clusters 
<tr><td> TkrNumHitsOnTracks         
<td>I<td>   Total number of TKR clusters on tracks
<tr><td> TkrFirstLayer
<td>I<td>   First layer containing a cluster 
<tr><td> TkrNumLayersHit
<td>I<td>   Total number of hit layers 
<tr><td> TkrNumGhosts         
<td>I<td>   Total number of TKR ghost clusters 
<tr><td> TkrNumDiags        
<td>I<td>   Total number of TKR ghost clusters as identified by TEM diagnostic info 
<tr><td> TkrNumBoth        
<td>I<td>   Total number of TKR clusters as regular and diagnostic ghosts 
<tr><td> TkrNumToT255s         
<td>I<td>   Total number of clusters with ToT==255
<tr><td> TkrNumSaturated
<td>I<td>   Number of clusters marked as saturated (all hits in each half-plane share the same ToT!)
<tr><td> TkrNumWideClusters         
<td>I<td>   Number of clusters 5 or more strips wide
<tr><td> TkrNumWiderClusters        
<td>I<td>   Number of clusters 9 or more strips wide
<tr><td> TkrNumSaturatedGhosts
<td>I<td>   Number of saturated ghost clusters
<tr><td> TkrNumWideGhosts
<td>I<td>   Number of wide ghost clusters
<tr><td> TkrNumWiderClusters        
<td>I<td>   Number of wider ghost clusters
<tr><td> TkrNumGhostsOnTracks         
<td>I<td>   Total number of ghost clusters on tracks
<tr><td> TkrNumDiagsOnTracks
<td>I<td>   Total number of track hits flagged as diagnostic ghosts
<tr><td> TkrNumBothOnTracks
<td>I<td>   Total number of hits flagged as both regular and diagnostic ghosts
<tr><td> TkrNumToT255sOnTracks         
<td>I<td>   Total number of ToT==255 clusters on tracks
<tr><td> TkrNumFlaggedTrackHits         
<td>I<td>   Total number of track hits flagged because
the track had some 255 or ghost hits
<tr><td> TkrNumSaturatedOnTracks
<td>I<td>   Number of saturated clusters on tracks
<tr><td> TkrNumWideClustersOnTracks        
<td>I<td>   Number of wide clusters on tracks
<tr><td> TkrNumWiderClustersOnTracks        
<td>I<td>   Number of wider clusters on tracks
<tr><td> TkrNumSaturatedGhostsOnTracks
<td>I<td>   Number of ghost clusters with saturated ToT
<tr><td> TkrNumWideGhostsOnTracks
<td>I<td>   Number of wide ghost clusters on tracks
<tr><td> TkrNumWiderGhostsOnTracks
<td>I<td>   Number of wider ghost clusters on tracks
<tr><td> TkrHitsInLyrNN, NN=(00,17)   
<td>I<td>   Number of clusters in (bi)layer NN 
(numbered from the bottom of the tracker) 
<tr><td> TkrNumRCTruncated
<td>I<td> Number of planes truncated due to RC buffers
<tr><td> TkrNumCCTruncated
<td>I<td> Number of planes truncated due to CC buffers

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

    addItem("TkrNumHits",             &Tkr_Cnv_Lyr_Hits);       
    addItem("TkrNumHitsOnTracks",     &Tkr_numHitsOnTracks);
    addItem("TkrFirstLayer",          &Tkr_Fst_Cnv_Lyr);        
    addItem("TkrNumLayersHit",        &Tkr_NCnv_Lyrs_Hit);
    addItem("TkrNumGhosts",           &Tkr_numGhosts);
    addItem("TkrNumDiags",            &Tkr_numDiags);
    addItem("TkrNumBoth",             &Tkr_numBoth);
    addItem("TkrNumToT255s",          &Tkr_numToT255s);
    addItem("TkrNumSaturated",        &Tkr_numSaturated);
    addItem("TkrNumWideClusters",     &Tkr_numWideClusters);
    addItem("TkrNumWiderClusters",    &Tkr_numWiderClusters);
    addItem("TkrNumSaturatedGhosts",        &Tkr_numSaturatedGhosts);
    addItem("TkrNumWideGhosts",      &Tkr_numWideGhostClusters);
    addItem("TkrNumWiderGhosts",     &Tkr_numWiderGhostClusters);
    addItem("TkrNumGhostsOnTracks",    &Tkr_numGhostsOnTracks);
    addItem("TkrNumDiagsOnTracks",     &Tkr_numDiagsOnTracks);
    addItem("TkrNumBothOnTracks",     &Tkr_numBothOnTracks);
    addItem("TkrNumToT255sOnTracks",    &Tkr_numToT255sOnTracks);
    addItem("TkrNumFlaggedTrackHits",   &Tkr_numFlaggedTrackHits);
    addItem("TkrNumSaturatedOnTracks", &Tkr_numSaturatedOnTracks);
    addItem("TkrNumWideClustersOnTracks",  &Tkr_numWideClustersOnTracks);
    addItem("TkrNumWiderClustersOnTracks", &Tkr_numWiderClustersOnTracks);
    addItem("TkrNumSaturatedGhostsOnTracks", &Tkr_numSaturatedGhostsOnTracks);
    addItem("TkrNumWideGhostsOnTracks",  &Tkr_numWideGhostClustersOnTracks);
    addItem("TkrNumWiderGhostsOnTracks", &Tkr_numWiderGhostClustersOnTracks);

    int i;
    char buffer[20];
    for(i=0;i<_nLayers;++i) {
        // vars of form TkrHitsInLyrNN, NN = (00,17), for _nLayers==18
        sprintf(buffer, "TkrHitsInLyr%02i",i);
        addItem(buffer, &Tkr_HitsPerLyr[i]);
    }

    // count truncation records
    addItem("TkrNumRCTruncs", &Tkr_numRCTruncs);
    addItem("TkrNumCCTruncs", &Tkr_numCCTruncs);

    zeroVals();

    return sc;
}


StatusCode TkrHitValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    //Make sure we have valid cluster data

    if (!pClusters) return sc;

    int layerIdx;
    for(layerIdx=0;layerIdx<_nLayers;++layerIdx) {
        int hitCount = 
            m_clusTool->getClusters(idents::TkrId::eMeasureX,layerIdx).size()
            + m_clusTool->getClusters(idents::TkrId::eMeasureY,layerIdx).size();

        if (hitCount > 0)
        {
            Tkr_Fst_Cnv_Lyr    = layerIdx;
            Tkr_NCnv_Lyrs_Hit += 1;
        }

        Tkr_HitsPerLyr[layerIdx] = hitCount;
    }

    Tkr_Cnv_Lyr_Hits = pClusters->size();

    Event::TkrClusterColConItr iter = pClusters->begin();
    for(; iter!=pClusters->end();++iter) {
        bool isGhost, isDiagGhost, isBoth, is255, isSaturated, isWide, isWider, isMarked;
        Event::TkrCluster* clust = *iter;

        is255       = clust->isSet(Event::TkrCluster::mask255);
        isGhost     = clust->isSet(Event::TkrCluster::maskGHOST);
        isDiagGhost = clust->isSet(Event::TkrCluster::maskDIAGNOSTIC);
        isBoth      = isGhost&&isDiagGhost;
        isSaturated = (clust->getRawToT()==250);
        isWide      = (clust->size()>=m_minWide);
        isWider     = (clust->size()>=m_minWider);
        isMarked = 
            clust->isSet(Event::TkrCluster::maskZAPGHOSTS);

        if(isMarked) {
            if     (is255)       Tkr_numToT255s++;
            else if(isBoth)      Tkr_numBoth++;
            else if(isGhost)     Tkr_numGhosts++;
            else if(isDiagGhost) Tkr_numDiags++;
            if(isSaturated) Tkr_numSaturatedGhosts++;
            if(isWide)      Tkr_numWideGhostClusters++;
            if(isWider)     Tkr_numWiderGhostClusters++;
        } else {
            if(isSaturated) Tkr_numSaturated++;
            if(isWide)      Tkr_numWideClusters++;
            if(isWider)     Tkr_numWiderClusters++;
        }

        bool onTrack = clust->hitFlagged();
        if(onTrack) {
            Tkr_numHitsOnTracks++;
            if(isMarked) {
                if(is255)            Tkr_numToT255sOnTracks++;
                else if(isBoth)      Tkr_numBothOnTracks++;
                else if(isGhost)     Tkr_numGhostsOnTracks++;
                else if(isDiagGhost) Tkr_numDiagsOnTracks++;
                else                 Tkr_numFlaggedTrackHits++;
                if(isSaturated) Tkr_numSaturatedGhostsOnTracks++;
                if(isWide)      Tkr_numWideGhostClustersOnTracks++;
                if(isWider)     Tkr_numWiderGhostClustersOnTracks++;
            }
            else {
                if(isSaturated) Tkr_numSaturatedOnTracks++;
                if(isWide)      Tkr_numWideClustersOnTracks++;
                if(isWider)     Tkr_numWiderClustersOnTracks++;
            }
        }
    }

    // truncation info
    SmartDataPtr<Event::TkrTruncationInfo>   
        pTruncInfo(m_pEventSvc,EventModel::TkrRecon::TkrTruncationInfo);
    //Make sure we have valid cluster data

    if (pTruncInfo) {
        Tkr_numRCTruncs = pTruncInfo->getNumRCTruncated();
        Tkr_numCCTruncs = pTruncInfo->getNumCCTruncated();
    }

    return sc;
}
