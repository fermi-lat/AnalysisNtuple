/** @file McTkrHitValsTool.cxx
    @brief declartion, implementaion of the class UserAlg

    $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/McTkrHitValsTool.cxx,v 1.6 2007/04/06 15:19:55 usher Exp $
*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McRelTableDefs.h"

#include "Event/Digi/TkrDigi.h"

// Volume identifier stuff
#include "idents/VolumeIdentifier.h"

#include <algorithm>
#include <set>

/*! @class McTkrHitValsTool
@brief calculates Monte Carlo values

@authors Bill Atwood, Leon Rochester
*/

class McTkrHitValsTool : public ValBase
{
public:
    
    McTkrHitValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~McTkrHitValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    // Internal functions
    int  GetTrackHits(const Event::McParticle* paricle);
    void CntTotalHits(const Event::McParticle* primary);
    int  GetSharedHits(const Event::McParticle* daught1, const Event::McParticle* daught2);
    void CntMcPosHits(const Event::McParticle* primary);

    // Variables for primary track
    int          m_primType;          // Particle type of primary
    int          m_primNumHits;       // Primary number of tracker hits 
                                      // (will be zero for a gamma)

    // Variables for first daughter
    int          m_dght1Type;         // Particle type of "best" daughter of primary
    int          m_dght1NumHits;      // "Best" daughter number of track hits

    // Variables for second daughter
    int          m_dght2Type;         // Particle type of "next" daughter track
    int          m_dght2NumHits;      // "Next" daughter number of track hits

    // McPositionHit/particle counts
    int          m_posHitPrimary;     // Number of McPositionHits associated to primary
    int          m_posHitGammas;      // Number of McPositionHits from gammas
    int          m_posHitElectrons;   // Number of McPositionHits not primary due to electrons
    int          m_posHitPositrons;   // Number of McPositionHits not primary due to positrons
    int          m_posHitOthers;      // Number of McPositionHits none of the above

    // Other variables
    int          m_process;           // Decay process for primary (if one)
    int          m_totalHits;         // Total number of non-noise hits in tracker
    int          m_totalPrimHits;     // Number of hits due to primary and its descendants
    int          m_totalPrimClusHits; // Number of clusters due to primary and its descendants
    int          m_dghtSharedHits;    // Number of clusters shared by the two daughters above

    // to decode the particle charge
    IParticlePropertySvc*  m_ppsvc;

    // Relational tables to be set up each event...
    Event::McPartToTrajectoryTab* m_partToTrajTab;
    Event::McPointToPosHitTab*    m_pointToPosHitTab;
    //Event::McPointToIntHitTab*    m_pointToIntHitTab;  // Future expansion
    Event::ClusMcPosHitTab*       m_clusToPosHitTab;
};

// Static factory for instantiation of algtool objects
static ToolFactory<McTkrHitValsTool> s_factory;
const IToolFactory& McTkrHitValsToolFactory = s_factory;

// Standard Constructor
McTkrHitValsTool::McTkrHitValsTool(const std::string& type, 
                               const std::string& name, 
                               const IInterface* parent)
                               : ValBase( type, name, parent ),
                                 m_partToTrajTab(0),
                                 m_pointToPosHitTab(0),
                                 m_clusToPosHitTab(0)
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}


/** @page anatup_vars_optional 
@section McTkrHitValstool McTkrHitValsTool Variables

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> McTHPrimType <td> I <td> Primary particle type (id)
<tr><td> McTHPrimNumHits <td> I <td> Primary number of Tracker hits (making clusters) 
<tr><td> McTHDght1Type <td> I <td> Daughter 1 (most hits) particle type (id)
<tr><td> McTHDght1NumHits <td> I <td> Daughter 1 number of Tracker hits
<tr><td> McTHDght2Type <td> I <td> Daughter 2 (most hits) particle type (id)
<tr><td> McTHDght2NumHits <td> I <td> Daughter 2 number of Tracker hits
<tr><td> mcTHPosHitPrimary <td> I <td> Number of McPositionHits associated to primary 
<tr><td> mcTHPosHitGamma <td> I <td> Total Number of McPositionHits due to gammas
<tr><td> mcTHPosHitElectron <td> I <td> Number of McPositionHits due to electrons
<tr><td> mcTHPosHitPositron <td> I <td> Total Number of McPositionHits due to positrons
<tr><td> mcTHPosHitOthers <td> I <td> Number of McPositionHits none of the above
<tr><td> mcTHTotalHits <td> I <td> Total number of MC generated Tracker hits 
</table>

*/


StatusCode McTkrHitValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;
  
    // get the services    
    if( serviceLocator() ) {
        if( service("ParticlePropertySvc", m_ppsvc, true).isFailure() ) {
            log << MSG::ERROR << "Service [ParticlePropertySvc] not found" << endreq;
        }
    } else {
        return StatusCode::FAILURE;
    }
 
	addItem("McTHPrimType",       &m_primType);
    addItem("McTHProcess",        &m_process);
	addItem("McTHPrimNumHits",    &m_primNumHits);
	addItem("McTHDght1Type",      &m_dght1Type);
	addItem("McTHDght1NumHits",   &m_dght1NumHits);
	addItem("McTHDght2Type",      &m_dght2Type);
	addItem("McTHDght2NumHits",   &m_dght2NumHits);
    addItem("McTHPosHitPrimary",  &m_posHitPrimary);
    addItem("McTHPosHitGamma",    &m_posHitGammas);
    addItem("McTHPosHitElectron", &m_posHitElectrons);
    addItem("McTHPosHitPositron", &m_posHitPositrons);
    addItem("McTHPosHitOthers",   &m_posHitOthers);
    addItem("McTHTotalHits",      &m_totalHits);
    addItem("McTHTotPrmHits",     &m_totalPrimHits);
    addItem("McTHTotPrmClusHits", &m_totalPrimClusHits);
    addItem("McTHDghtSharedHits", &m_dghtSharedHits);

    zeroVals();
    
    return sc;
}

StatusCode McTkrHitValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream  log( msgSvc(), name() );

    // Grab the event header to get the event time for printing a message in debug mode
    SmartDataPtr<Event::EventHeader>   header(m_pEventSvc,    EventModel::EventHeader);
    double t = header->time();
    log << MSG::DEBUG << "Event time: " << t << endreq;;

    // First task is to recover the McParticle Collection and categorize the event
    SmartDataPtr<Event::McParticleCol> particleCol(m_pEventSvc, EventModel::MC::McParticleCol);

    // No collection then we do not have MC! 
    if (particleCol == 0) return sc;

    // Recover Relational tables that we will need...
    SmartDataPtr<Event::McPartToTrajectoryTabList> 
        mcPartToTraj(m_pEventSvc, "/Event/MC/McPartToTrajectory");
    m_partToTrajTab = new Event::McPartToTrajectoryTab(mcPartToTraj);

    // Recover the McPositionHits to trajectory points relational tables
    SmartDataPtr<Event::McPointToPosHitTabList> 
        mcPointToPosHitList(m_pEventSvc, "/Event/MC/McPointToPosHit");
    m_pointToPosHitTab = new Event::McPointToPosHitTab(mcPointToPosHitList);

    // Recover the McIntegratingHits to trajectory points relational tables
    //SmartDataPtr<Event::McPointToIntHitTabList> 
    //    mcPointToIntHitList(m_pEventSvc, "/Event/MC/McPointToIntHit");
    //m_pointToIntHitTab = new Event::McPointToIntHitTab(mcPointToIntHitList);

    // Retrieve the table giving us TkrDigi to McPositionHit relations
    SmartDataPtr<Event::ClusMcPosHitTabList> 
        clusMcPosHitList(m_pEventSvc, EventModel::Digi::TkrClusterHitTab);

    // It can happen that there are no clusters...
    if (clusMcPosHitList != 0) m_clusToPosHitTab = new Event::ClusMcPosHitTab(clusMcPosHitList);

    // Our first task is to identify the primary in the event
    Event::McParticle* primary    = 0;
    std::string        primString = "primary";

    Event::McParticleCol::iterator mcPartIter = particleCol->begin();

    for(; mcPartIter != particleCol->end(); mcPartIter++)
    {
        Event::McParticle* mcPart = *mcPartIter;

        // Is it the primary?
        if (mcPart->getProcess() == primString) 
        {
            primary = mcPart;
            break;
        }
    }

    // Everything follows if we have found a primary
    if (primary)
    {
        // We now count the total number of hits in the Tracker and the number
        // due to the primary and its products
        CntTotalHits(primary);

        // Count the number of McPositionHits due to primary, electrons and positrons
        CntMcPosHits(primary);

        // Next up is to get the number of hits for each of McParticles in the list
        // Set up the map to hold the results
        std::map<const Event::McParticle*,int> mcPartHitMap;
        mcPartHitMap.clear();

        // Reset the McParticle iterator
        mcPartIter = particleCol->begin();

        for(++mcPartIter; mcPartIter != particleCol->end(); mcPartIter++)
        {
            Event::McParticle* mcPart = *mcPartIter;

            // Get the number of hits from this track 
            int numHits = GetTrackHits(mcPart);

            mcPartHitMap[mcPart] = numHits;
        }

        // Recover the primary's particle identification and number of hits
        m_primType    = primary->particleProperty();
        m_primNumHits = mcPartHitMap[primary];

        // Get some of the primary's characteristics
//        ParticleProperty* ppty     = m_ppsvc->findByStdHepID(m_primType);
//        std::string       partName = ppty->particle(); 

        // We will need pointers to the "two" daughters to determine shared hits
        const Event::McParticle* daughter1 = 0;
        const Event::McParticle* daughter2 = 0;

        // Loop through the daughters keeping track of the number of hits
        const SmartRefVector<Event::McParticle>& daughterList = primary->daughterList();

        for(SmartRefVector<Event::McParticle>::const_iterator daughterIter = daughterList.begin();
            daughterIter != daughterList.end(); daughterIter++)
        {
            const Event::McParticle* daughter = *daughterIter;

            int numHits = mcPartHitMap[daughter];

            // Check to see if this exceeds the current daughter1 count
            if (numHits > m_dght1NumHits || m_dght1Type == 0)
            {
                // Check to see if we need to switch down to second daughter
                if (m_dght1NumHits > m_dght2NumHits || m_dght2Type == 0)
                {
                    m_dght2NumHits = m_dght1NumHits;
                    m_dght2Type    = m_dght2Type;
                    daughter2      = daughter1;
                }

                m_dght1NumHits = numHits;
                m_dght1Type    = daughter->particleProperty();
                daughter1      = daughter;
            }
            // Otherwise, check to see if this is the hit count for the next particle
            else if (numHits > m_dght2NumHits || m_dght2Type == 0)
            {
                m_dght2NumHits = numHits;
                m_dght2Type    = daughter->particleProperty();
                daughter2      = daughter;
            }
        }

        // Finally... get the number of hits these two daughters share
        m_dghtSharedHits = GetSharedHits(daughter1, daughter2);
    }

    // Clean up after ourselves
    if (m_partToTrajTab)    {delete m_partToTrajTab;    m_partToTrajTab    = 0;}
    if (m_pointToPosHitTab) {delete m_pointToPosHitTab; m_pointToPosHitTab = 0;}
    //if (m_pointToIntHitTab) {delete m_pointToIntHitTab; m_pointToIntHitTab = 0;}
    if (m_clusToPosHitTab)  {delete m_clusToPosHitTab;  m_clusToPosHitTab  = 0;}
    
    log << MSG::DEBUG << " returning. " << endreq;

    return sc;
}


void McTkrHitValsTool::CntTotalHits(const Event::McParticle* primary)
{
    m_totalHits     = 0;
    m_totalPrimHits = 0;

    // If there is no cluster table then there are no hits (by definition)
    if (m_clusToPosHitTab)
    {
        // Retrieve pointer to McPositionHit Collection in the TDS
        SmartDataPtr<Event::McPositionHitCol> posHitCol(m_pEventSvc, EventModel::MC::McPositionHitCol);

        // Define the results containers
        std::map<int, std::set<Event::TkrCluster*> > idToHitMap;    // For each track, number of actual hits made
        std::set<int>                                primTrkIdSet;  // For primary, number of tracks resulting
        std::set<Event::TkrCluster*>                 primClusSet;   // Number of clusters by tracks from primary

        idToHitMap.clear();
        primTrkIdSet.clear();
        primClusSet.clear();

        // Loop through the McPositionHit collection
        for (Event::McPositionHitCol::iterator posHitIter = posHitCol->begin(); posHitIter != posHitCol->end(); posHitIter++)
        {
            Event::McPositionHit* posHit = *posHitIter;

            // Is this hit related to a Cluster (ie did this make a hit in the Tracker?)
            Event::ClusMcPosHitVec clusHitVec = m_clusToPosHitTab->getRelBySecond(posHit);

            // Do we have anything? (note that this will eliminate McPositionHits from the ACD)
           
            if (!clusHitVec.empty())
            {
                // An McPositionHit belongs to one cluster (but a cluster can be made up of
                // several McPositionHits...)
                Event::TkrCluster* cluster = (*(clusHitVec.begin()))->getFirst();

                // Check our current map for this hit, start by getting the track id for this hit
                int trackId = posHit->getPackedFlags(); // Track ID is hiding in here...

                idToHitMap[trackId].insert(cluster);

                // Does this particle descend from the primary?
                //if (posHit->originMcParticle() == primary) // ok, this turns out to not be valid for gammas... 
                if (posHit->originMcParticle() != 0)         // temporary kludge...
                {
                    primClusSet.insert(cluster);
                    primTrkIdSet.insert(trackId);
                }
            }
        }

        // Now pass through our map to count total number of hits from each track 
        for (std::map<int,std::set<Event::TkrCluster*> >::iterator idToHitIter = idToHitMap.begin(); 
             idToHitIter != idToHitMap.end(); idToHitIter++)
        {
            m_totalHits += idToHitIter->second.size();
        }

        // Go through and count total number of hits from the primary and its descendants
        for (std::set<int>::iterator trkIdIter = primTrkIdSet.begin(); trkIdIter != primTrkIdSet.end(); trkIdIter++)
        {
            m_totalPrimHits += idToHitMap[*trkIdIter].size();
        }

        // And now the total number of clusters from the primary and its descendants
        m_totalPrimClusHits = primClusSet.size();
    }

    return;
}

int McTkrHitValsTool::GetTrackHits(const Event::McParticle* particle)
{
    int nTrackHits = 0;

    // If there is no cluster table then there are no hits
    if (m_clusToPosHitTab)
    {
        // Look up trajectory for this particle
        Event::McPartToTrajectoryVec trajVec = m_partToTrajTab->getRelByFirst(particle);
        Event::McTrajectory* mcTraj          = (*(trajVec.begin()))->getSecond();

        // Set up to loop through Trajectory points
        std::vector<Event::McTrajectoryPoint*> points = mcTraj->getPoints();
        std::vector<Event::McTrajectoryPoint*>::const_iterator pointIter;

        // Use a set to keep track of unique clusters (since a cluster may be associated
        // to several McPositionHits). We'll use this to count Tracker hits
        std::set<Event::TkrCluster*> trackClusters;
        trackClusters.clear();

        // Loop through all the McTrajectoryPoints on this track
        for(pointIter = points.begin(); pointIter != points.end(); pointIter++) 
        {
            // De-reference the McTrajectoryPoint pointer
            Event::McTrajectoryPoint* mcPoint = *pointIter;

            // Look for McTrajectoryPoint <--> McPositionHit
            Event::McPointToPosHitVec posRelVec = m_pointToPosHitTab->getRelByFirst(mcPoint);

            // Did this trajectory point leave a mark?
            if (!posRelVec.empty())
            {
                // Can only be one McPositionHit per McTrajectoryPoint... by definition
                Event::McPointToPosHitRel* hitRel = *(posRelVec.begin());

                Event::McPositionHit* mcPosHit = hitRel->getSecond();

                // Get any TkrClusters associated to this McPositionHit (note: this eliminates ACD)
                Event::ClusMcPosHitVec clusHitVec = m_clusToPosHitTab->getRelBySecond(mcPosHit);

                // Do we have anything?
                if (!clusHitVec.empty())
                {
                    // An McPositionHit belongs to one cluster (but a cluster can be made up of
                    // several McPositionHits...)
                    Event::TkrCluster* cluster = (*(clusHitVec.begin()))->getFirst();

                    // "cluster" will be unique at the end 
                    trackClusters.insert(cluster);
                }
            }
        }

        // ok, the size of trackClusters is, then, the number of hits this track left
        nTrackHits = trackClusters.size();
    }

    return nTrackHits;
}

int McTkrHitValsTool::GetSharedHits(const Event::McParticle* daughter1, const Event::McParticle* daughter2)
{
    int numShared = 0;

    // Plan? 
    // loop through hits on daughter1 track and get clusters. 
    //      use table to get mcposhits from cluster
    //      loop through mcposhits and look for match to daughter2

    // If there is no cluster table then there are no hits, also protect that both daughters exist
    if (m_clusToPosHitTab && daughter1 && daughter2)
    {
        // Look up trajectory for this particle 
        // (use the second daughter because it is shorter)
        Event::McPartToTrajectoryVec trajVec = m_partToTrajTab->getRelByFirst(daughter2);
        Event::McTrajectory* mcTraj          = (*(trajVec.begin()))->getSecond();

        // Set up to loop through Trajectory points
        std::vector<Event::McTrajectoryPoint*> points = mcTraj->getPoints();
        std::vector<Event::McTrajectoryPoint*>::const_iterator pointIter;

        // We need a list of clusters resulting from the daughter track, we will then try to relate
        // those clusters to those used by the other daughter track 
        std::set<Event::TkrCluster*> trackClusters;
        trackClusters.clear();

        // Loop through all the McTrajectoryPoints on this track and store away the clusters
        for(pointIter = points.begin(); pointIter != points.end(); pointIter++) 
        {
            // De-reference the McTrajectoryPoint pointer
            Event::McTrajectoryPoint* mcPoint = *pointIter;

            // Look for McTrajectoryPoint <--> McPositionHit
            Event::McPointToPosHitVec posRelVec = m_pointToPosHitTab->getRelByFirst(mcPoint);

            // Did this trajectory point leave a mark?
            if (!posRelVec.empty())
            {
                // Can only be one McPositionHit per McTrajectoryPoint... by definition
                Event::McPointToPosHitRel* hitRel = *(posRelVec.begin());

                Event::McPositionHit* mcPosHit = hitRel->getSecond();

                // Get any TkrClusters associated to this McPositionHit (note: this eliminates ACD)
                Event::ClusMcPosHitVec clusHitVec = m_clusToPosHitTab->getRelBySecond(mcPosHit);

                // Do we have anything?
                if (!clusHitVec.empty())
                {
                    // An McPositionHit belongs to one cluster (but a cluster can be made up of
                    // several McPositionHits...)
                    Event::TkrCluster* cluster = (*(clusHitVec.begin()))->getFirst();

                    // "cluster" will be unique at the end 
                    trackClusters.insert(cluster);
                }
            }
        }

        // Go through the set of clusters to look up any relation they might have to the other track
        for (std::set<Event::TkrCluster*>::iterator clusIter = trackClusters.begin(); clusIter != trackClusters.end(); clusIter++)
        {
            Event::TkrCluster* cluster = *clusIter;

            // Look up all McPositionHits associated with this cluster
            Event::ClusMcPosHitVec clusHitVec = m_clusToPosHitTab->getRelByFirst(cluster);

            // Proceed if it might be shared
            if (clusHitVec.size() > 1)
            {
                for(Event::ClusMcPosHitVec::iterator relHitIter = clusHitVec.begin(); relHitIter != clusHitVec.end(); relHitIter++)
                {
                    Event::McPositionHit* relPosHit = (*relHitIter)->getSecond();

                    // Check that a valid pointer exists and see if it matches the first daughter
                    if (relPosHit && relPosHit->mcParticle() == daughter1)
                    {
                        numShared++;
                        break;
                    }
                }
            }
        }
    }

    return numShared;
}

void McTkrHitValsTool::CntMcPosHits(const Event::McParticle* primary)
{
    // Define codes for electrons and positrons
    static Event::McParticle::StdHepId gamma    =  22;
    static Event::McParticle::StdHepId electron =  11;
    static Event::McParticle::StdHepId positron = -11;

    // Clear counters
    m_posHitPrimary   = 0;
    m_posHitGammas    = 0;
    m_posHitElectrons = 0;
    m_posHitPositrons = 0;
    m_posHitOthers    = 0;

    // First task is to recover the McParticle Collection and categorize the event
    SmartDataPtr<Event::McPositionHitCol> positionHitCol(m_pEventSvc, EventModel::MC::McPositionHitCol);

    // If collection then proceed
    if (positionHitCol)
    {
        // Loop over the collection of McPositionHits
        for(Event::McPositionHitCol::iterator posHitIter = positionHitCol->begin();
            posHitIter != positionHitCol->end(); posHitIter++)
        {
            Event::McPositionHit* posHit = *posHitIter;

            // Check the volume identifier to cut on tracker hits
            idents::VolumeIdentifier volId = posHit->volumeID();

            // Do not count ACD hits
            if (volId[0] != 0) continue;

            // Get parent information
            const Event::McParticle* particle = posHit->mcParticle();
            const Event::McParticle* mother   = 0;
            
            if (particle) mother = particle->getMother();

            // Check to see if hit is associated to primary or immediate daughters
            if      (particle                  == primary)  m_posHitPrimary++;
            else if (mother                    == primary)  m_posHitPrimary++;

            // What kind of particle made the hit?
            if      (posHit->getMcParticleId() == gamma)    m_posHitGammas++;
            else if (posHit->getMcParticleId() == electron) m_posHitElectrons++;
            else if (posHit->getMcParticleId() == positron) m_posHitPositrons++;
            else                                            m_posHitOthers++;
        }
    }


    return;
}
