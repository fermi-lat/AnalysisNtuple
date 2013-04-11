/** @file TreeValsTool.cxx
@brief Calculates the Tracker Tree variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/TreeValsTool.cxx,v 1.21 2013/01/18 22:38:26 usher Exp $
*/

// Include files
#include "AnalysisNtuple/IValsTool.h"
#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"

//*** Special version to look at tree based tracking timing
#include "GaudiKernel/IChronoStatSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrFilterParams.h"
#include "Event/Recon/TreeClusterRelation.h"
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
    int getNumLeavesThisBranch(const Event::TkrVecNode* node, int curNumLeaves);

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

    IPropagator*           m_G4PropTool; 

    IChronoStatSvc*        m_chronoSvc;

    // Tree stuff (if present)
    float Tkr_num_vecPoints;
    float Tkr_ave_vecPoints;
    float Tkr_lyr_wVecPoints;
    float Tkr_lyr_mostVecPoints;
    float Tkr_num_mostVecPoints;

    float Tkr_num_vecLinks;
    float Tkr_num_trees;
    float Tkr_tree1_nTrks;
    float Tkr_tree1_depth;
    float Tkr_tree1_nLeaves;

    float Tkr_tree1_nNodes;
        float Tkr_tree1_Thin_nNodes;
        float Tkr_tree1_Thick_nNodes;
        float Tkr_tree1_Thin_Leaves;
        float Tkr_tree1_Thick_Leaves;
        float Tkr_tree1_Blank_nNodes;
        float Tkr_tree1_ThinNodes_RLn;
        float Tkr_tree1_ThickNodes_RLn;

    float Tkr_tree1_bestBranchAngToAxis;
    float Tkr_tree1_axisSeededAngToAxis;

    float Tkr_tree1_nBranches;
    float Tkr_tree1_firstLinkAngle;
    float Tkr_tree1_firstLinkSum;
    float Tkr_tree1_bestRms;
    float Tkr_tree1_bestBiLayers;
    float Tkr_tree1_maxWidthLyr;
    float Tkr_tree1_maxWidth;
    float Tkr_tree1_lastWidth;
    float Tkr_tree1_bestLeaves;
    float Tkr_tree1_scndLeaves;
    float Tkr_tree1_scndRms;
    float Tkr_tree1_scndDepth;
    float Tkr_tree1_PosX;
    float Tkr_tree1_PosY;
    float Tkr_tree1_PosZ;
    float Tkr_tree1_DirX;
    float Tkr_tree1_DirY;
    float Tkr_tree1_DirZ;
    float Tkr_tree1_NumBoxes;
    float Tkr_tree1_ChiSquare;
    float Tkr_tree1_RmsTrans;
    float Tkr_tree1_RmsLong;

    float Tkr_tree1_numLinks;
    float Tkr_tree1_calPosDoca;
    float Tkr_tree1_calDoca68;
    float Tkr_tree1_calDoca95;
    float Tkr_tree1_calDocaMax;

    float Tkr_tree2_nTrks;
    float Tkr_tree2_depth;
    float Tkr_tree2_nLeaves;
    float Tkr_tree2_nNodes;
    float Tkr_tree2_nBranches;
    float Tkr_tree2_bestRms;
    float Tkr_tree2_bestBiLayers;
    float Tkr_tree2_maxWidthLyr;
    float Tkr_tree2_maxWidth;
    float Tkr_tree2_lastWidth;

    // For now, hide the TkrFilterParams output down here... 
    float TFP_numParams;
    float TFP_bestPosX;
    float TFP_bestPosY;
    float TFP_bestPosZ;
    float TFP_bestDirX;
    float TFP_bestDirY;
    float TFP_bestDirZ;
    float TFP_bestNumHits;
    float TFP_bestChiSquare;
    float TFP_bestAveDist;
    float TFP_bestRmsTrans;
    float TFP_bestRmsLong;
    float TFP_bestRmsLongAsym;

    float TFP_numLinks;
    float TFP_calDoca68;
    float TFP_calDoca95;
    float TFP_calDocaMax;

    // For now, the timing info
    float Aud_treeBasedTool;
    float Aud_treeBasedTool_link;
    float Aud_treeBasedTool_node;
    float Aud_treeBasedTool_build;
    float Aud_tkrVecLinkBuilderTool_singleLink;
    float Aud_tkrVecLinkBuilderTool_multiLink;
    float Aud_houghFilterTool;
    float Aud_houghFilterTool_fill;
    float Aud_houghFilterTool_peak;
    float Aud_houghFilterTool_build;
};

// Static factory for instantiation of algtool objects
//static ToolFactory<TreeValsTool> s_factory;
//const IToolFactory& TreeValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TreeValsTool);

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

        if (service("ChronoStatSvc", m_chronoSvc, true).isFailure())
        {
            log << MSG::ERROR << "Couldn't find the ChronoSvc!" << endreq;
            return StatusCode::FAILURE;
        }

    } else {
        return fail;
    }

    // load up the map

    // Tree stuff
    addItem("TkrNumVecPoints",             &Tkr_num_vecPoints);
    addItem("TkrAveVecPtsLyr",             &Tkr_ave_vecPoints);
    addItem("TkrLyrsWVecPoints",           &Tkr_lyr_wVecPoints);
    addItem("TkrLyrMostVecPoints",         &Tkr_lyr_mostVecPoints);
    addItem("TkrNumMostVecPoints",         &Tkr_num_mostVecPoints);
                                     
    addItem("TkrNumVecLinks",              &Tkr_num_vecLinks);
    addItem("TkrNumTrees",                 &Tkr_num_trees);
    addItem("TkrTree1NTrks",               &Tkr_tree1_nTrks);
    addItem("TkrTree1Depth",               &Tkr_tree1_depth);
    addItem("TkrTree1Leaves",              &Tkr_tree1_nLeaves);
                                     
    addItem("TkrTree1Nodes",               &Tkr_tree1_nNodes);
        addItem("TkrTree1ThinNodes",           &Tkr_tree1_Thin_nNodes);
        addItem("TkrTree1ThickNodes",          &Tkr_tree1_Thick_nNodes);
        addItem("TkrTree1ThinLeaves",          &Tkr_tree1_Thin_Leaves);
        addItem("TkrTree1ThickLeaves",         &Tkr_tree1_Thick_Leaves);
                                     
    addItem("TkrTree1BlankNodes",          &Tkr_tree1_Blank_nNodes);
        addItem("TkrTree1ThinRLnNodes",        &Tkr_tree1_ThinNodes_RLn);
        addItem("TkrTree1ThickRLnNodes",       &Tkr_tree1_ThickNodes_RLn);

    addItem("TkrTree1BestBranchAngToAxis", &Tkr_tree1_bestBranchAngToAxis);
    addItem("TkrTree1AxisSeededAngToAxis", &Tkr_tree1_axisSeededAngToAxis);
                                     
    addItem("TkrTree1Branches",            &Tkr_tree1_nBranches);
    addItem("TkrTree1FirstLinkAngle",      &Tkr_tree1_firstLinkAngle);
    addItem("TkrTree1FirstLinkSum",        &Tkr_tree1_firstLinkSum);
    addItem("TkrTree1BestRms",             &Tkr_tree1_bestRms);
    addItem("TkrTree1BestBiLayers",        &Tkr_tree1_bestBiLayers);
    addItem("TkrTree1MaxWidthLyr",         &Tkr_tree1_maxWidthLyr);
    addItem("TkrTree1MaxWidth",            &Tkr_tree1_maxWidth);
    addItem("TkrTree1LastWidth",           &Tkr_tree1_lastWidth);
    addItem("TkrTree1BestLeaves",          &Tkr_tree1_bestLeaves);
    addItem("TkrTree1ScndLeaves",          &Tkr_tree1_scndLeaves);
    addItem("TkrTree1ScndRms",             &Tkr_tree1_scndRms);
    addItem("TkrTree1ScndDepth",           &Tkr_tree1_scndDepth);
    addItem("TkrTree1PosX",                &Tkr_tree1_PosX);
    addItem("TkrTree1PosY",                &Tkr_tree1_PosY);
    addItem("TkrTree1PosZ",                &Tkr_tree1_PosZ);
    addItem("TkrTree1DirX",                &Tkr_tree1_DirX);
    addItem("TkrTree1DirY",                &Tkr_tree1_DirY);
    addItem("TkrTree1DirZ",                &Tkr_tree1_DirZ);
    addItem("TkrTree1NumBoxes",            &Tkr_tree1_NumBoxes);
    addItem("TkrTree1ChiSquare",           &Tkr_tree1_ChiSquare);
    addItem("TkrTree1RmsTrans",            &Tkr_tree1_RmsTrans);
    addItem("TkrTree1RmsLong",             &Tkr_tree1_RmsLong);

    addItem("TkrTree1NumLinks",            &Tkr_tree1_numLinks);
    addItem("TkrTree1CalPosDoca",          &Tkr_tree1_calPosDoca);
    addItem("TkrTree1CalDoca68",           &Tkr_tree1_calDoca68);
    addItem("TkrTree1CalDoca95",           &Tkr_tree1_calDoca95);
    addItem("TkrTree1CalDocaMax",          &Tkr_tree1_calDocaMax);
                                     
    addItem("TkrTree2NTrks",               &Tkr_tree2_nTrks);
    addItem("TkrTree2Depth",               &Tkr_tree2_depth);
    addItem("TkrTree2Leaves",              &Tkr_tree2_nLeaves);
    addItem("TkrTree2Nodes",               &Tkr_tree2_nNodes);
    addItem("TkrTree2Branches",            &Tkr_tree2_nBranches);
    addItem("TkrTree2BestRms",             &Tkr_tree2_bestRms);
    addItem("TkrTree2BestBiLayers",        &Tkr_tree2_bestBiLayers);
    addItem("TkrTree2MaxWidthLyr",         &Tkr_tree2_maxWidthLyr);
    addItem("TkrTree2MaxWidth",            &Tkr_tree2_maxWidth);
    addItem("TkrTree2LastWidth",           &Tkr_tree2_lastWidth);
                                     
    addItem("TFPNumParams",                &TFP_numParams);
    addItem("TFPBestPosX",                 &TFP_bestPosX);
    addItem("TFPBestPosY",                 &TFP_bestPosY);
    addItem("TFPBestPosZ",                 &TFP_bestPosZ);
    addItem("TFPBestDirX",                 &TFP_bestDirX);
    addItem("TFPBestDirY",                 &TFP_bestDirY);
    addItem("TFPBestDirZ",                 &TFP_bestDirZ);
    addItem("TFPBestNumHits",              &TFP_bestNumHits);
    addItem("TFPBestChiSquare",            &TFP_bestChiSquare);
    addItem("TFPBestAverageDist",          &TFP_bestAveDist);
    addItem("TFPBestRmsTrans",             &TFP_bestRmsTrans);
    addItem("TFPBestRmsLong",              &TFP_bestRmsLong);
    addItem("TFPBestRmsLongAsym",          &TFP_bestRmsLongAsym);

    addItem("TFPNumLinks",                 &TFP_numLinks);
    addItem("TFPCalDoca68",                &TFP_calDoca68);
    addItem("TFPCalDoca95",                &TFP_calDoca95);
    addItem("TFPCalDocaMax",               &TFP_calDocaMax);
                                     
    addItem("AudTreeBasedTool",            &Aud_treeBasedTool);
    addItem("AudTreeBasedTool_link",       &Aud_treeBasedTool_link);
    addItem("AudTreeBasedTool_node",       &Aud_treeBasedTool_node);
    addItem("AudTreeBasedTool_build",      &Aud_treeBasedTool_build);
    addItem("AudLinkTool_singleLink",      &Aud_tkrVecLinkBuilderTool_singleLink);
    addItem("AudLinkTool_multiLink",       &Aud_tkrVecLinkBuilderTool_multiLink);
    addItem("AudHoughFilterTool",          &Aud_houghFilterTool);
    addItem("AudHoughFilterTool_fill",     &Aud_houghFilterTool_fill);
    addItem("AudHoughFilterTool_peak",     &Aud_houghFilterTool_peak);
    addItem("AudHoughFilterTool_build",    &Aud_houghFilterTool_build);

    zeroVals();

    return sc;
}

StatusCode TreeValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // How many TkrVecPoints?
    Event::TkrVecPointInfo* vecPointInfo = 
        SmartDataPtr<Event::TkrVecPointInfo>(m_pEventSvc, EventModel::TkrRecon::TkrVecPointInfo);

    if (vecPointInfo) 
    {
        Tkr_num_vecPoints = vecPointInfo->getNumTkrVecPoints();

        // Go through and determine the average number of vecpoints per layer
        int nLyrsWVecPoints =  0;
        int lyrWMostPoints  = -1;
        int nLyrWMostPoints =  0;

        for(Event::TkrLyrToVecPointItrMapItr vecMapItr  = vecPointInfo->getLyrToVecPointItrMap()->begin();
                                             vecMapItr != vecPointInfo->getLyrToVecPointItrMap()->end();
                                             vecMapItr++)
        {
            int numLinks = std::distance(vecMapItr->second.first, vecMapItr->second.second);

            if (numLinks) 
            {
                nLyrsWVecPoints++;

                if (numLinks > nLyrWMostPoints)
                {
                    nLyrWMostPoints = numLinks;
                    lyrWMostPoints  = vecMapItr->first;
                }
            }
        }

        Tkr_lyr_wVecPoints = nLyrsWVecPoints;

        if (nLyrsWVecPoints > 0)
        {
            Tkr_ave_vecPoints     = Tkr_num_vecPoints / Tkr_lyr_wVecPoints;
            Tkr_lyr_mostVecPoints = lyrWMostPoints;
            Tkr_num_mostVecPoints = nLyrWMostPoints;
        }
    }

    // How many TkrVecPointLinks?
    SmartDataPtr<Event::TkrVecPointsLinkCol> vecLinksCol(m_pEventSvc,EventModel::TkrRecon::TkrVecPointsLinkCol);

    if (vecLinksCol) Tkr_num_vecLinks = vecLinksCol->size();

    // Do tree stuff here
    SmartDataPtr<Event::TkrTreeCol> treeCol(m_pEventSvc,"/Event/TkrRecon/TkrTreeCol");

    if (treeCol)
    {
        // Recover the map between tree and relations to calorimeter
        Event::TreeToRelationMap* treeToRelationMap = SmartDataPtr<Event::TreeToRelationMap>(m_pEventSvc, EventModel::Recon::TreeToRelationMap);
        Event::CalCluster*        calCluster        = 0;

        Tkr_num_trees = treeCol->size();

        if (Tkr_num_trees > 0 && (*treeCol->begin())->getHeadNode())
        {
            Event::TkrTreeColConPtr  treeItr  = treeCol->begin(); 
            const Event::TkrTree*    tree     = *treeItr++;
            const Event::TkrVecNode* headNode = tree->getHeadNode();
            const Event::TkrVecNode* bestLeaf = tree->getBestLeaf();
            const Event::TkrVecNode* scndLeaf = tree->getSecondLeaf();

            Tkr_tree1_nTrks                   = tree->size();
            Tkr_tree1_depth                   = headNode->getDepth();
            Tkr_tree1_nLeaves                 = headNode->getNumLeaves();

            Tkr_tree1_nNodes                  = headNode->getNumNodesInTree();
                        Tkr_tree1_Thin_nNodes             = headNode->getNumThinNodesInTree();
                        Tkr_tree1_Thick_nNodes            = headNode->getNumThickNodesInTree();
                        Tkr_tree1_Thin_Leaves             = headNode->getNumThinLeavesInTree();
                        Tkr_tree1_Thick_Leaves            = headNode->getNumThickLeavesInTree();
                        Tkr_tree1_Blank_nNodes            = headNode->getNumBlankNodesInTree();
                        Tkr_tree1_ThinNodes_RLn           = headNode->getNumThinRLnInTree();
                        Tkr_tree1_ThickNodes_RLn          = headNode->getNumThickRLnInTree();

            Tkr_tree1_bestBranchAngToAxis     = tree->getBestBranchAngleToAxis();
            Tkr_tree1_axisSeededAngToAxis     = tree->getAxisSeededAngleToAxis();
                        
            Tkr_tree1_nBranches               = headNode->getNumBranches();
            Tkr_tree1_bestRms                 = headNode->getBestRmsAngle();
            Tkr_tree1_bestBiLayers            = headNode->getBestNumBiLayers();
            Tkr_tree1_firstLinkAngle          = 0.;

            if (!headNode->front()->empty())
            {
                const Event::TkrVecNode* firstNode = headNode->front();
                const Event::TkrVecNode* scndNode  = firstNode->front();

                double cosLinkAngle = 
                    firstNode->getAssociatedLink()->getVector().dot(scndNode->getAssociatedLink()->getVector());

                cosLinkAngle             = std::max(std::min(cosLinkAngle, 1.), -1.);
                Tkr_tree1_firstLinkAngle = acos(cosLinkAngle);
                Tkr_tree1_firstLinkSum   = sqrt(scndNode->getRmsAngleSum());
            }

            if (bestLeaf) Tkr_tree1_bestLeaves   = getNumLeavesThisBranch(bestLeaf, 0);

            if (scndLeaf)
            {
                Tkr_tree1_scndLeaves   = getNumLeavesThisBranch(scndLeaf, 0);
                Tkr_tree1_scndDepth    = scndLeaf->getTreeStartLayer() - scndLeaf->getCurrentBiLayer() + 1;
                Tkr_tree1_scndRms      = scndLeaf->getBestRmsAngle();
            }

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

            // Try calculating link doca to cal cluster
            // To do that we need to use the relations between tree and cluster, provided it exists
            if (treeToRelationMap)
            {
                Event::TreeToRelationMap::iterator relationItr = treeToRelationMap->find(const_cast<Event::TkrTree*>(tree));

                if (relationItr != treeToRelationMap->end())
                {

                    calCluster = relationItr->second.front()->getCluster();
                }

                // If an associated cluster then calculate the docas
                if (calCluster)
                {
                  // Ph.Bruel: switch to centroid postion corrected forthe hodoscopic effect
                  //                  const Point& clusPos = calCluster->getPosition();
                  Point clusPos = calCluster->getCorPosition(tree->getAxisParams()->getEventAxis());

                    // Now we need to visit EVERY link in the tree. To do so we use a list...
                    std::list<const Event::TkrVecNode*> nodeList;

                    // Add the first node to the list to start to populate it
                    nodeList.push_back(tree->getHeadNode());

                    // Create an iterator to run over the entries in the list
                    std::list<const Event::TkrVecNode*>::const_iterator nodeListItr = nodeList.begin();

                    // Create a vector to keep track of the docas we calculate
                    std::vector<double> treeDocaVec;
					std::vector<double> filterDocaVec;

                    // We loop until we exhaust the entries in the list
                    while(nodeListItr != nodeList.end())
                    {
                        const Event::TkrVecNode* node = *nodeListItr;

                        // First lets add the child nodes of this node to our list
                        for(Event::TkrVecNodeSet::const_iterator childItr = node->begin(); childItr != node->end(); childItr++)
                        {
                            nodeList.push_back(*childItr);
                        }

                        // Ok, now increment the node list iterator
                        nodeListItr++;

                        // Recover the associated link (if there is one)
                        const Event::TkrVecPointsLink* link = node->getAssociatedLink();

                        if (link)
                        {
                            // We want to determine the distance of the link, projected to a disk perpendicular to the 
                            // axis from the link position to the cluster centroid, at the cluster centroid. 
                            // Start with the axis from the link position to the cluster centroid
                            Vector linkToPos = clusPos - link->getPosition();

                            // Get the angle between this axis and the link
                            double cosTheta  = link->getVector().dot(linkToPos.unit());

                            // Compute the arclength along the link direction to the disk
                            double arcLen    = linkToPos.magnitude() / cosTheta;

                            // Calculate the point on this disk
                            Point  docaPos   = link->getPosition() + arcLen * link->getVector();

                            // Get a vector from the cluster centroid to this point
                            Vector clus2Doca = docaPos - clusPos;

                            // Get the doca
                            double doca      = clus2Doca.magnitude();

                            treeDocaVec.push_back(doca);

							// Was this link used by the filter?
							if (link->getStatusBits() & 0x03000000) filterDocaVec.push_back(doca);
                        }
                    }

                    // Sort out vector of docas
                    std::sort(treeDocaVec.begin(), treeDocaVec.end());

                    // Extract the 68%, 95% and final elements
                    int idx68 = int(0.68 * float(treeDocaVec.size()));
                    int idx95 = int(0.95 * float(treeDocaVec.size()));

                    Tkr_tree1_numLinks   = treeDocaVec.size();
                    Tkr_tree1_calPosDoca = relationItr->second.front()->getTreeClusDoca();
                    Tkr_tree1_calDoca68  = treeDocaVec[idx68];
                    Tkr_tree1_calDoca95  = treeDocaVec[idx95];
                    Tkr_tree1_calDocaMax = treeDocaVec.back();

					// Since we are here, check the docas for the filter as well. 
					// But only if there is something there
					if (!filterDocaVec.empty())
					{

					// Since we are here, check the docas for the filter as well. 
					std::sort(filterDocaVec.begin(), filterDocaVec.end());

                    // Extract the 68%, 95% and final elements
                    idx68 = int(0.68 * float(filterDocaVec.size()));
                    idx95 = int(0.95 * float(filterDocaVec.size()));

					// Fill the TFP variables
					TFP_numLinks   = filterDocaVec.size();
					TFP_calDoca68  = filterDocaVec[idx68];
					TFP_calDoca95  = filterDocaVec[idx95];
					TFP_calDocaMax = filterDocaVec.back();
					}
                }
            }

            Tkr_tree1_maxWidthLyr = firstLayer - Tkr_tree1_maxWidthLyr;

            // Get the axis parameters
            const Event::TkrFilterParams* treeAxis = tree->getAxisParams();

            if (treeAxis)
            {
                const Point filterPos  = treeAxis->getEventPosition();
                const Vector filterDir = treeAxis->getEventAxis();

                Tkr_tree1_PosX         = filterPos.x();
                Tkr_tree1_PosY         = filterPos.y();
                Tkr_tree1_PosZ         = filterPos.z();
                                   
                Tkr_tree1_DirX         = filterDir.x();
                Tkr_tree1_DirY         = filterDir.y();
                Tkr_tree1_DirZ         = filterDir.z();
                                   
                Tkr_tree1_NumBoxes     = treeAxis->getNumHitsTotal();
                Tkr_tree1_ChiSquare    = treeAxis->getChiSquare();
                Tkr_tree1_RmsTrans     = treeAxis->getTransRms();
                Tkr_tree1_RmsLong      = treeAxis->getLongRms();
            }

            if (treeItr != treeCol->end())
            {
                tree     = *treeItr++;
                headNode = tree->getHeadNode();

                Tkr_tree2_nTrks        = tree->size();
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

    // Do tree stuff here
    SmartDataPtr<Event::TkrFilterParamsCol> filterParamsCol(m_pEventSvc,"/Event/TkrRecon/TkrFilterParamsCol");

    if (filterParamsCol)
    {
        TFP_numParams = filterParamsCol->size();

        if (TFP_numParams > 0)
        {
            Event::TkrFilterParams* filterParams = filterParamsCol->front();

            Point filterPos     = filterParams->getEventPosition();
            Vector filterDir    = filterParams->getEventAxis();

            TFP_bestPosX        = filterPos.x();
            TFP_bestPosY        = filterPos.y();
            TFP_bestPosZ        = filterPos.z();

            TFP_bestDirX        = filterDir.x();
            TFP_bestDirY        = filterDir.y();
            TFP_bestDirZ        = filterDir.z();

            TFP_bestNumHits     = filterParams->getNumHitsTotal();
            TFP_bestChiSquare   = filterParams->getChiSquare();
            TFP_bestAveDist     = filterParams->getAverageDistance();
            TFP_bestRmsTrans    = filterParams->getTransRms();
            TFP_bestRmsLong     = filterParams->getLongRms();
            TFP_bestRmsLongAsym = filterParams->getLongRmsAsym();
        }
    }

    // Get the auditor info
    IChronoStatSvc::ChronoTime time = m_chronoSvc->chronoDelta("TreeBasedTool",IChronoStatSvc::USER);
    Aud_treeBasedTool = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TreeBasedTool_link",IChronoStatSvc::USER);
    Aud_treeBasedTool_link = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TreeBasedTool_node",IChronoStatSvc::USER);
    Aud_treeBasedTool_node = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TreeBasedTool_build",IChronoStatSvc::USER);
    Aud_treeBasedTool_build = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TkrVecLinkBuilderTool_singleLink",IChronoStatSvc::USER);
    Aud_tkrVecLinkBuilderTool_singleLink = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TkrVecLinkBuilderTool_multiLink",IChronoStatSvc::USER);
    Aud_tkrVecLinkBuilderTool_multiLink = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TkrHoughFilterTool",IChronoStatSvc::USER);
    Aud_houghFilterTool = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TkrHoughFilterTool_fill",IChronoStatSvc::USER);
    Aud_houghFilterTool_fill = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TkrHoughFilterTool_peak",IChronoStatSvc::USER);
    Aud_houghFilterTool_peak = static_cast<float>(time)*0.000001;
    time = m_chronoSvc->chronoDelta("TkrHoughFilterTool_build",IChronoStatSvc::USER);
    Aud_houghFilterTool_build = static_cast<float>(time)*0.000001;

    return sc;
}

int TreeValsTool::getNumLeavesThisBranch(const Event::TkrVecNode* node, int curNumLeaves)
{
    // We are tracing up from the leaf to the tree parent
    // If we are on either the "best" branch, or the "next" best branch, but not both,
    // then we can simply count the number of leaves to and below this node
    if (!(node->isOnBestBranch() && node->isOnNextBestBranch()))
    {
        curNumLeaves = node->getNumLeaves();
    }
    // Otherwise, we are on the branch that eventually splits into both branches... at this point we 
    // need to count the leaves emanating from this branch to add to the current count of curNumLeaves
    else
    {
        if (!node->empty())
        {
            Event::TkrVecNodeSet::const_iterator nodeItr = node->begin();

            // Loop over nodes but skip the first which is always going to be the "best"
            for(++nodeItr; nodeItr != node->end(); nodeItr++)
            {
                // Be sure to not count the number of leaves on the best or second branches
                if (!((*nodeItr)->isOnBestBranch() || (*nodeItr)->isOnNextBestBranch())) 
                    curNumLeaves += (*nodeItr)->getNumLeaves();
            }
        }
    }

    // If we have a valid parent then counting is not done yet! 

    if (node->getParentNode()) curNumLeaves = getNumLeavesThisBranch(node->getParentNode(), curNumLeaves);

    return curNumLeaves;
}
