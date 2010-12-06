/** @file TreeValsTool.cxx
@brief Calculates the Tracker Tree variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TreeValsTool.cxx,v 1.2 2010/11/24 20:01:56 usher Exp $
*/

// Include files
#include "AnalysisNtuple/IValsTool.h"
#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrTree.h"
#include "geometry/Ray.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrFlagHitsTool.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GaudiKernel/IToolSvc.h"
#include "Doca.h"

#include <cstring>


// M_PI defined in ValBase.h

/*! @class TreeValsTool
@brief calculates Tkr values

@authors Bill Atwood, Leon Rochester
*/

class TreeValsTool :  public ValBase 
{
public:

    TreeValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~TreeValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

private:

    double towerEdge(Point pos) const;
    double containedFraction(Point pos, double gap, double r, 
        double costh, double phi) const;
    float SSDEvaluation(const Event::TkrTrack* track); 

    // some local constants
    double m_towerPitch;
    int    m_xNum;
    int    m_yNum;
    double m_activeWidth;
    bool   m_useNew;
    bool   m_enableVetoDiagnostics;
    int    m_messageCount;

    // some pointers to services

    /// pointer to tracker geometry
    ITkrGeometrySvc*       m_tkrGeom;
    /// GlastDetSvc used for access to detector info
    IGlastDetSvc*          m_detSvc; 
    /// pointer to queryclusterstool
    ITkrQueryClustersTool* m_queryClusters;
    /// 
//    IPropagatorSvc* m_propSvc;

    IPropagator* m_G4PropTool; 

    // Tree stuff (if present)
    float Tkr_num_vecPoints;
    float Tkr_num_vecLinks;
    float Tkr_num_trees;
    float Tkr_tree1_depth;
    float Tkr_tree1_nLeaves;
    float Tkr_tree1_nNodes;
    float Tkr_tree1_nBranches;
    float Tkr_tree1_bestRms;
    float Tkr_tree1_bestBiLayers;
    float Tkr_tree1_maxWidthLyr;
    float Tkr_tree1_maxWidth;
    float Tkr_tree1_lastWidth;
    float Tkr_tree2_depth;
    float Tkr_tree2_nLeaves;
    float Tkr_tree2_nNodes;
    float Tkr_tree2_nBranches;
    float Tkr_tree2_bestRms;
    float Tkr_tree2_bestBiLayers;
    float Tkr_tree2_maxWidthLyr;
    float Tkr_tree2_maxWidth;
    float Tkr_tree2_lastWidth;
};

// Static factory for instantiation of algtool objects
static ToolFactory<TreeValsTool> s_factory;
const IToolFactory& TreeValsToolFactory = s_factory;

// Standard Constructor
TreeValsTool::TreeValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this);

    //declareProperty("useNew", m_useNew=true);
    // Old SSDVeto variable has been removed; always use "new" code
    m_useNew = true;
//    declareProperty("enableVetoDiagnostics", m_enableVetoDiagnostics=false);

}

//************************
// Below to be updated for tree stuff when it settles down
//************************

/** @page anatup_vars 
@section TreeValsTool TreeValsTool Variables

Notes: 
- Variables called Tkr1Xxx refer to the "best" track; 
those called Tkr2Xxx refer to the second track.
- A number of variables have the word "Hits" in their name. This <em>usually</em>
refers to clusters! In some cases it refers to TkrTrackHits. 
The description should make it clear which meaning is intended.
- For variables listed as Tkr[1/2]Xxx there are two versions in the ntuple, 
one for the best and one for the second track. 
- The labels are not entirely consistent, but it's probably 
too disruptive to fix them at this point.
For example: TkrRadLength, TkrTrackLength, TkrTwrEdge refer to track 1. 
Also, Tkr2Angle and Tkr2HDoca are quantities that depend on both tracks.
- The variables associated with the second track are undefined 
if there is only one track! 
Check TkrNumTracks before using these variables! 
In fact check TkrNumTracks before using first-track variables, 
for the same reason.
- A new section of (optional) ssd-veto diagnostic variables has been added. 
They are not written out by default.
- Several new variables starting with "TkrV" have been added. These refer to 
quantities associated with the track likely to have cause the Acd veto.
- Some deleted variables, all Tkr2: FirstHits, DifHits, Gaps, FirstGaps,
DieEdge, KalThetaMs, [X/Y/Z]Dir, Phi, Theta, [X/Y/Z]0.

@subsection general General variables
<table>
<tr><th> Variable <th> Type  <th> Description				                 
<tr><td> TkrNumTracks 	
<td>F<td>   Number of tracks found (Maximum is set by TkrRecon, currently 10) 
</table>

@subsection both Variables that exist for both best and second tracks

<table>
<tr><th> Variable <th> Type  <th> Description				                 
<tr><td> Tkr[1/2]Chisq 
<td>F<td>   Track chisquared 
<tr><td> Tkr[1/2]FirstChisq  
<td>F<td>   Track chisquared for first Tkr[1/2]FirstHits layers  
</table>

*/

StatusCode TreeValsTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;

    MsgStream log(msgSvc(), name());

    log << MSG::INFO  << "#################" << endreq << "# ";
    log << (m_useNew ? "New " : "Old ");
    log << "version" << endreq << "#################" << endreq;

    if((ValBase::initialize()).isFailure()) return StatusCode::FAILURE;

    m_messageCount = 0;

    // get the services

    if( serviceLocator()) {

        if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return fail;
        }

        m_towerPitch = m_tkrGeom->towerPitch();
        m_xNum       = m_tkrGeom->numXTowers();
        m_yNum       = m_tkrGeom->numYTowers();
        m_activeWidth = m_tkrGeom->nWaferAcross()*m_tkrGeom->waferPitch() + 
            (m_tkrGeom->nWaferAcross()-1)*m_tkrGeom->ladderGap();

        // find GlastDevSvc service
        if (service("GlastDetSvc", m_detSvc, true).isFailure()){
            log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
            return StatusCode::FAILURE;
        }

        IToolSvc* toolSvc = 0;
        if(service("ToolSvc", toolSvc, true).isFailure()) {
            log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
            return StatusCode::FAILURE;
        }
        if(!toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)) {
            log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
            return StatusCode::FAILURE;
        }

    } else {
        return fail;
    }

    // load up the map

    // Tree stuff
    addItem("TkrNumVecPoints",      &Tkr_num_vecPoints);
    addItem("TkrNumVecLinks",       &Tkr_num_vecLinks);
    addItem("TkrNumTrees",          &Tkr_num_trees);
    addItem("TkrTree1Depth",        &Tkr_tree1_depth);
    addItem("TkrTree1Leaves",       &Tkr_tree1_nLeaves);
    addItem("TkrTree1Nodes",        &Tkr_tree1_nNodes);
    addItem("TkrTree1Branches",     &Tkr_tree1_nBranches);
    addItem("TkrTree1BestRms",      &Tkr_tree1_bestRms);
    addItem("TkrTree1BestBiLayers", &Tkr_tree1_bestBiLayers);
    addItem("TkrTree1MaxWidthLyr",  &Tkr_tree1_maxWidthLyr);
    addItem("TkrTree1MaxWidth",     &Tkr_tree1_maxWidth);
    addItem("TkrTree1LastWidth",    &Tkr_tree1_lastWidth);
    addItem("TkrTree2Depth",        &Tkr_tree2_depth);
    addItem("TkrTree2Leaves",       &Tkr_tree2_nLeaves);
    addItem("TkrTree2Nodes",        &Tkr_tree2_nNodes);
    addItem("TkrTree2Branches",     &Tkr_tree2_nBranches);
    addItem("TkrTree2BestRms",      &Tkr_tree2_bestRms);
    addItem("TkrTree2BestBiLayers", &Tkr_tree2_bestBiLayers);
    addItem("TkrTree2MaxWidthLyr",  &Tkr_tree2_maxWidthLyr);
    addItem("TkrTree2MaxWidth",     &Tkr_tree2_maxWidth);
    addItem("TkrTree2LastWidth",    &Tkr_tree2_lastWidth);

    zeroVals();

    return sc;
}

StatusCode TreeValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // How many TkrVecPoints?
    SmartDataPtr<Event::TkrVecPointCol> vecPointCol(m_pEventSvc,EventModel::TkrRecon::TkrVecPointCol);

    if (vecPointCol) Tkr_num_vecPoints = vecPointCol->size();

    // How many TkrVecPointLinks?
    SmartDataPtr<Event::TkrVecPointsLinkCol> vecLinksCol(m_pEventSvc,EventModel::TkrRecon::TkrVecPointsLinkCol);

    if (vecLinksCol) Tkr_num_vecLinks = vecLinksCol->size();

    // Do tree stuff here
    SmartDataPtr<Event::TkrTreeCol> treeCol(m_pEventSvc,"/Event/TkrRecon/TkrTreeCol");

    if (treeCol)
    {
        Tkr_num_trees = treeCol->size();

        if (Tkr_num_trees > 0)
        {
            Event::TkrTreeColConPtr  treeItr  = treeCol->begin(); 
            const Event::TkrTree*    tree     = *treeItr++;
            const Event::TkrVecNode* headNode = tree->getHeadNode();

            Tkr_tree1_depth        = headNode->getDepth();
            Tkr_tree1_nLeaves      = headNode->getNumLeaves();
            Tkr_tree1_nNodes       = headNode->getNumNodesInTree();
            Tkr_tree1_nBranches    = headNode->getNumBranches();
            Tkr_tree1_bestRms      = headNode->getBestRmsAngle();
            Tkr_tree1_bestBiLayers = headNode->getBestNumBiLayers();

            Event::TkrNodeSiblingMap::const_reverse_iterator sibItr = tree->getSiblingMap()->rbegin();

            Tkr_tree1_maxWidthLyr  = sibItr->first;
            Tkr_tree1_maxWidth     = sibItr->second.size();
            Tkr_tree1_lastWidth    = tree->getSiblingMap()->begin()->second.size();

            
            float firstLayer = sibItr->first;
    
            for(; sibItr != tree->getSiblingMap()->rend(); sibItr++)
            {
                const std::vector<const Event::TkrVecNode*>& nodeVec = sibItr->second;

                int nodeVecSize = nodeVec.size();

                if (nodeVecSize > Tkr_tree1_maxWidth)
                {
                    Tkr_tree1_maxWidthLyr = sibItr->first;
                    Tkr_tree1_maxWidth    = nodeVecSize;
                }
            }

            Tkr_tree1_maxWidthLyr = firstLayer - Tkr_tree1_maxWidthLyr;

            if (treeItr != treeCol->end())
            {
                tree     = *treeItr++;
                headNode = tree->getHeadNode();

                Tkr_tree2_depth        = headNode->getDepth();
                Tkr_tree2_nLeaves      = headNode->getNumLeaves();
                Tkr_tree2_nNodes       = headNode->getNumNodesInTree();
                Tkr_tree2_nBranches    = headNode->getNumBranches();
                Tkr_tree2_bestRms      = headNode->getBestRmsAngle();
                Tkr_tree2_bestBiLayers = headNode->getBestNumBiLayers();

                sibItr = tree->getSiblingMap()->rbegin();

                Tkr_tree2_maxWidthLyr  = sibItr->first;
                Tkr_tree2_maxWidth     = sibItr->second.size();
                Tkr_tree2_lastWidth    = tree->getSiblingMap()->begin()->second.size();

                float firstLayer = sibItr->first;
    
                for(; sibItr != tree->getSiblingMap()->rend(); sibItr++)
                {
                    const std::vector<const Event::TkrVecNode*>& nodeVec = sibItr->second;

                    int nodeVecSize = nodeVec.size();

                    if (nodeVecSize > Tkr_tree2_maxWidth)
                    {
                        Tkr_tree2_maxWidthLyr = sibItr->first;
                        Tkr_tree2_maxWidth    = nodeVecSize;
                    }
                }

                Tkr_tree2_maxWidthLyr = firstLayer - Tkr_tree2_maxWidthLyr;
            }
        }
    }

    return sc;
}
