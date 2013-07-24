/** @file McTkrHitValsTool.cxx
    @brief declartion, implementaion of the class UserAlg

    $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/McTkrHitValsTool.cxx,v 1.19 2013/07/04 17:34:39 carmelo Exp $
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
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McRelTableDefs.h"

#include "Event/Digi/TkrDigi.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"

// Volume identifier stuff
#include "idents/VolumeIdentifier.h"
#include "idents/AcdId.h"
#include "idents/CalXtalId.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include <algorithm>
#include <set>

#include "Doca.h"

//cs test prop
#include "GaudiKernel/IToolSvc.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"

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
    int  GetAcdHits(const Event::McParticle* paricle);
    int  GetTrackHits(const Event::McParticle* paricle);
    int  GetCalHits(const Event::McParticle* paricle);
    void CntTotalHits(const Event::McParticle* primary);
    int  GetSharedHits(const Event::McParticle* daught1, const Event::McParticle* daught2);
    void CntMcPosHits(const Event::McParticle* primary);
    int  GetInteracted();

    //Attempt to count the number of truncated planes in the MC particle trajectory
    void countTruncatedPlanes(Event::McParticle* primary);

    // Variables for primary track
    int          m_primType;          // Particle type of primary
    int          m_numDaughters;      // Number of daughters for the primary
    int          m_primNumTkrHits;    // Primary number of tracker hits 
    int          m_primNumAcdHits;    // Primary number of Acd hits 
    int          m_primNumCalHits;    // Primary number of calorimeter hits 
    int          m_primInteracted;    // Primary interacted in LAT
    int          m_primInteractFlag;             // Primary first interaction flag: 0 for conversion, 1 compton, 2 other
    int          m_primInteractMaterialId;       // Primary first interaction material Id.
    char         m_primInteractMaterialName[80]; // Primary first interaction material Name (string)
    float        m_primInteractX0;    // Primary first interaction X0
    float        m_primInteractY0;    // Primary first interaction Y0
    float        m_primInteractZ0;    // Primary first interaction Z0
    float        m_primInteractEleXDir;
    float        m_primInteractEleYDir;
    float        m_primInteractEleZDir;
    float        m_primInteractEleEne;
    float        m_primInteractPosXDir;
    float        m_primInteractPosYDir;
    float        m_primInteractPosZDir;
    float        m_primInteractPosEne;

    // Variables for first daughter
    int          m_dght1Type;         // Particle type of "best" daughter of primary
    int          m_dght1NumTkrHits;   // "Best" daughter number of track hits
    int          m_dght1NumAcdHits;   // "Best" daughter number of acd hits
    int          m_dght1NumCalHits;   // "Best" daughter number of cal hits
    char         m_dght1Process[80];  // Decay process for this daughter

    // Variables for second daughter
    int          m_dght2Type;         // Particle type of "next" daughter track
    int          m_dght2NumTkrHits;   // "Next" daughter number of track hits
    int          m_dght2NumAcdHits;   // "Next" daughter number of acd hits
    int          m_dght2NumCalHits;   // "Next" daughter number of cal hits
    char         m_dght2Process[80];  // Decay process for this daughter

    // McPositionHit/particle counts
    int          m_posHitPrimary;     // Number of McPositionHits associated to primary
    int          m_posHitGammas;      // Number of McPositionHits from gammas
    int          m_posHitElectrons;   // Number of McPositionHits not primary due to electrons
    int          m_posHitPositrons;   // Number of McPositionHits not primary due to positrons
    int          m_posHitOthers;      // Number of McPositionHits none of the above

    // Other variables
    int          m_totalTkrHits;      // Total number of non-noise hits in tracker
    int          m_totalPrimTkrHits;  // Number of tracker hits due to primary and its descendants
    int          m_totalPrimAcdHits;  // Number of acd hits due to primary and its descendants
    int          m_totalPrimCalHits;  // Number of cal hits due to primary and its descendants
    int          m_totalPrimClusHits; // Number of clusters due to primary and its descendants
    int          m_dghtSharedHits;    // Number of clusters shared by the two daughters above
    float        m_totEnergyDepCal;   // Total energy deposited in the calorimeter
    float        m_totEnergyDepAcd;   // Total energy deposited in the ACD

    // Variables for keeping track of truncations
    int          m_planesCrossed;     // Number of planes traversed by primary daughters
    int          m_topPlane;          // Top plane in event
    int          m_bottomPlane;       // Bottom plane in event
    int          m_nTruncatedPlanes;  // Number of truncated planes encountered
    int          m_nTruncatedLost;    // Number of planes with lost hits due to truncation
    float        m_truncActDist[36];  // Array of truncated plane "active distances"

    /// pointer to tracker geometry
    ITkrGeometrySvc*      m_tkrGeom;

    // to decode the particle charge
    IParticlePropertySvc* m_ppsvc;    

  IPropagator* m_G4PropTool;  //cs test prop

    // Useful tracker geometry constants
    float m_siThickness;
    float m_stripLength;
    float m_activeWidth;
    float m_towerPitch;
    int   m_numXTowers;
    int   m_numYTowers;

    // Relational tables to be set up each event...
    Event::McPartToTrajectoryTab*      m_partToTrajTab;
    Event::McPointToPosHitTab*         m_pointToPosHitTab;
    Event::McPointToIntHitTab*         m_pointToIntHitTab;  
    Event::ClusMcPosHitTab*            m_clusToPosHitTab;

    std::map<idents::CalXtalId, float> m_xtalIdToEnergyMap;
    std::map<idents::AcdId,     float> m_acdIdToEnergyMap;
};

// Static factory for instantiation of algtool objects
//static ToolFactory<McTkrHitValsTool> s_factory;
//const IToolFactory& McTkrHitValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(McTkrHitValsTool);

// Standard Constructor
McTkrHitValsTool::McTkrHitValsTool(const std::string& type, 
                               const std::string& name, 
                               const IInterface* parent)
                               : ValBase( type, name, parent ),
                                 m_partToTrajTab(0),
                                 m_pointToPosHitTab(0),
                                 m_pointToIntHitTab(0),
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
<tr><td> McTHPrimNumTkrHits <td> I <td> Primary number of Tracker hits (making clusters) 
<tr><td> McTHPrimNumAcdHits <td> I <td> Primary number of Acd hits (McPositionHits)
<tr><td> McTHPrimNumCalHits <td> I <td> Primary number of Cal hits (McIntegratingHits)
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
<tr><td> McTHPrimInteraction <td> I <td> Primary particle first interaction flag: 0 for conversion, 1 compton, 2 other
<tr><td> McTHPrimInterMatId <td> I <td> Id of the material in which the first interaction occurs
<tr><td> McTHPrimInterMatName <td> C <td> Name of the material in which the first interaction occurs
<tr><td> McTHPrimInter[X/Y/Z]0 <td> F <td> Position of the first interaction
<tr><td> McTHPrimInter[Ele/Pos][X/Y/Z]Dir <td> F <td> Initial direction of Electron/Positron pair
<tr><td> McTHPrimInter[Ele/Pos]E <td> F <td> Initial Energy of Electron/Positron pair
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
        if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return StatusCode::FAILURE;
        }
        
        //cs test prop
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
        return StatusCode::FAILURE;
    }
 
    addItem("McTHPrimType",         &m_primType);     
    addItem("McTHNumDaughters",     &m_numDaughters); 
    addItem("McTHPrimNumHits",      &m_primNumTkrHits);
    addItem("McTHPrimNumAcdHits",   &m_primNumAcdHits);
    addItem("McTHPrimNumCalHits",   &m_primNumCalHits);
    addItem("McTHPrimInteracted",   &m_primInteracted);
    addItem("McTHDght1Type",        &m_dght1Type);
    addItem("McTHDght1NumTkrHits",  &m_dght1NumTkrHits);
    addItem("McTHDght1NumAcdHits",  &m_dght1NumAcdHits);
    addItem("McTHDght1NumCalHits",  &m_dght1NumCalHits);
    addItem("McTHDght1Process",      m_dght1Process);
    addItem("McTHDght2Type",        &m_dght2Type);
    addItem("McTHDght2NumTkrHits",  &m_dght2NumTkrHits);
    addItem("McTHDght2NumAcdHits",  &m_dght2NumAcdHits);
    addItem("McTHDght2NumCalHits",  &m_dght2NumCalHits);
    addItem("McTHDght2Process",      m_dght2Process);
    addItem("McTHPosHitPrimary",    &m_posHitPrimary);
    addItem("McTHPosHitGamma",      &m_posHitGammas);
    addItem("McTHPosHitElectron",   &m_posHitElectrons);
    addItem("McTHPosHitPositron",   &m_posHitPositrons);
    addItem("McTHPosHitOthers",     &m_posHitOthers);
    addItem("McTHTotalTkrHits",     &m_totalTkrHits);
    addItem("McTHTotPrmTkrHits",    &m_totalPrimTkrHits);
    addItem("McTHTotPrmAcdHits",    &m_totalPrimAcdHits);
    addItem("McTHTotPrmCalHits",    &m_totalPrimCalHits);
    addItem("McTHTotPrmClusHits",   &m_totalPrimClusHits);
    addItem("McTHDghtSharedHits",   &m_dghtSharedHits);
    addItem("McTHTotEneDepCal",     &m_totEnergyDepCal);
    addItem("McTHTotEneDepAcd",     &m_totEnergyDepAcd);

    addItem("McTHPlanesCrossed",    &m_planesCrossed);
    addItem("McTHTopPlane",         &m_topPlane);
    addItem("McTHBottomPlane",      &m_bottomPlane);
    addItem("McTHnTruncatedPlanes", &m_nTruncatedPlanes);
    addItem("McTHnTruncatedLost",   &m_nTruncatedLost);

    addItem("McTHPrimInteraction",  &m_primInteractFlag);
    addItem("McTHPrimInterMatId",   &m_primInteractMaterialId);
    addItem("McTHPrimInterMatName",  m_primInteractMaterialName);
    addItem("McTHPrimInterX0",      &m_primInteractX0); // These 3 may duplicate other
    addItem("McTHPrimInterY0",      &m_primInteractY0); // quantities and may be
    addItem("McTHPrimInterZ0",      &m_primInteractZ0); // removed in the future
    addItem("McTHPrimInterEleXDir", &m_primInteractEleXDir);
    addItem("McTHPrimInterEleYDir", &m_primInteractEleYDir);
    addItem("McTHPrimInterEleZDir", &m_primInteractEleZDir);
    addItem("McTHPrimInterEleE",    &m_primInteractEleEne);
    addItem("McTHPrimInterPosXDir", &m_primInteractPosXDir);
    addItem("McTHPrimInterPosYDir", &m_primInteractPosYDir);
    addItem("McTHPrimInterPosZDir", &m_primInteractPosZDir);
    addItem("McTHPrimInterPosE",    &m_primInteractPosEne);

    char buffer[20];
    for(int idx = 0; idx < 36; idx++) 
    {
        // vars of form McTHTruncActDist_, NN = (00,35), for 36 planes
        sprintf(buffer, "McTHTruncActDist_%02i",idx);
        addItem(buffer, &m_truncActDist[idx]);
    }


    zeroVals();

    // Recover useful constants
    m_towerPitch  = m_tkrGeom->towerPitch();
    m_siThickness = m_tkrGeom->siThickness();
    m_numXTowers  = m_tkrGeom->numXTowers();
    m_numYTowers  = m_tkrGeom->numYTowers();
    m_stripLength = m_tkrGeom->trayWidth();
    m_activeWidth = m_tkrGeom->nWaferAcross()*m_tkrGeom->waferPitch() 
                  + (m_tkrGeom->nWaferAcross()-1)*m_tkrGeom->ladderGap();
    
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
    SmartDataPtr<Event::McPointToIntHitTabList> 
        mcPointToIntHitList(m_pEventSvc, "/Event/MC/McPointToIntHit");
    m_pointToIntHitTab = new Event::McPointToIntHitTab(mcPointToIntHitList);

    // Retrieve the table giving us TkrDigi to McPositionHit relations
    SmartDataPtr<Event::ClusMcPosHitTabList> 
        clusMcPosHitList(m_pEventSvc, EventModel::Digi::TkrClusterHitTab);

    // Clear energy deposited maps
    m_xtalIdToEnergyMap.clear();
    m_acdIdToEnergyMap.clear();

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
        // Primary's number of daughters
        m_numDaughters = primary->daughterList().size();

        // We now count the total number of hits in the Tracker and the number
        // due to the primary and its products
        CntTotalHits(primary);

        // Count the number of McPositionHits due to primary, electrons and positrons
        CntMcPosHits(primary);

        // At this point count the number of truncated planes
        countTruncatedPlanes(primary);

        // Next up is to get the number of hits in the various detectors for each of McParticles in the list
        // First set up maps to hold these hits for each detector element
        // Start with Tracker
        std::map<const Event::McParticle*,int> mcPartToTkrHitMap;
        mcPartToTkrHitMap.clear();

        // Then ACD 
        std::map<const Event::McParticle*,int> mcPartToAcdHitMap;
        mcPartToAcdHitMap.clear();

        // Then the CAL
        std::map<const Event::McParticle*,int> mcPartToCalHitMap;
        mcPartToCalHitMap.clear();

        // Reset the McParticle iterator
        mcPartIter = particleCol->begin();

        // Loop through the list of particles and build the variouls hit maps
        for(++mcPartIter; mcPartIter != particleCol->end(); mcPartIter++)
        {
            Event::McParticle* mcPart = *mcPartIter;

            // Get the number of TKR hits from this track 
            int numHits = GetTrackHits(mcPart);
            mcPartToTkrHitMap[mcPart] = numHits;

            // Get the number of ACD hits from this track 
            int numAcdHits = GetAcdHits(mcPart);
            mcPartToAcdHitMap[mcPart] = numAcdHits;

            // Get the number of CAL hits from this track 
            int numCalHits = GetCalHits(mcPart);
            mcPartToCalHitMap[mcPart] = numCalHits;
        }

        // Recover the primary's particle identification and number of hits
        m_primType         = primary->particleProperty();
        m_primNumTkrHits   = mcPartToTkrHitMap[primary];
        m_primNumAcdHits   = mcPartToAcdHitMap[primary];
        m_primNumCalHits   = mcPartToCalHitMap[primary];

        // We will need pointers to the "two" daughters to determine shared hits
        const Event::McParticle* daughter1 = 0;
        const Event::McParticle* daughter2 = 0;

        // Loop through the daughters keeping track of the number of hits
        const SmartRefVector<Event::McParticle>& daughterList = primary->daughterList();

        for(SmartRefVector<Event::McParticle>::const_iterator daughterIter = daughterList.begin();
            daughterIter != daughterList.end(); daughterIter++)
        {
            const Event::McParticle* daughter = *daughterIter;

            // Tracker hits first
            int numHits = mcPartToTkrHitMap[daughter];

            // Check to see if this exceeds the current daughter1 count
            if (numHits > m_dght1NumTkrHits || m_dght1Type == 0)
            {
                // Check to see if we need to switch down to second daughter
                if (m_dght1NumTkrHits > m_dght2NumTkrHits || m_dght2Type == 0)
                {
                    m_dght2NumTkrHits = m_dght1NumTkrHits;
                    m_dght2Type       = m_dght2Type;
                    daughter2         = daughter1;
                }

                m_dght1NumTkrHits = numHits;
                m_dght1Type       = daughter->particleProperty();
                daughter1         = daughter;
            }
            // Otherwise, check to see if this is the hit count for the next particle
            else if (numHits > m_dght2NumTkrHits || m_dght2Type == 0)
            {
                m_dght2NumTkrHits = numHits;
                m_dght2Type       = daughter->particleProperty();
                daughter2         = daughter;
            }
        }

        // ACD and CAL hits associate to daughter 1
        m_dght1NumAcdHits = mcPartToAcdHitMap[daughter1];
        m_dght1NumCalHits = mcPartToCalHitMap[daughter1];

        // ACD and CAL hits associate to daughter 2
        m_dght2NumAcdHits = mcPartToAcdHitMap[daughter2];
        m_dght2NumCalHits = mcPartToCalHitMap[daughter2];

        // Count the total ACD hits due to the simulation
        m_totalPrimAcdHits = m_acdIdToEnergyMap.size();

        // Total energy deposited in ACD
        for(std::map<idents::AcdId, float>::iterator acdMapIter  = m_acdIdToEnergyMap.begin();
                                                     acdMapIter != m_acdIdToEnergyMap.end();
                                                     acdMapIter++)
        {
            m_totEnergyDepAcd += acdMapIter->second;
        }

        // Now count the total CAL hits due to the simulation
        m_totalPrimCalHits = m_xtalIdToEnergyMap.size();

        // Total energy deposited in CAL
        for(std::map<idents::CalXtalId, float>::iterator calMapIter  = m_xtalIdToEnergyMap.begin();
                                                         calMapIter != m_xtalIdToEnergyMap.end();
                                                         calMapIter++)
        {
            m_totEnergyDepCal += calMapIter->second;
        }

        // Recover the decay process for this primary (which is conveniently located in its daughter)
        std::string process1 = daughter1 ? daughter1->getProcess() : "";
        strncpy(m_dght1Process, process1.c_str(),80);
        if (m_dght1Process == "") memset(m_dght1Process, '_', 80);

        std::string process2 = daughter2 ? daughter2->getProcess() : "";
        strncpy(m_dght2Process, process2.c_str(),80);
        if (m_dght2Process == "") memset(m_dght2Process, '_', 80);

        // Finally... get the number of hits these two daughters share
        m_dghtSharedHits = GetSharedHits(daughter1, daughter2);


        // Now some polarization-specific quantities for photons: look for the first interaction of the photon.

        // init output vars
        m_primInteractMaterialId = -1; 
        std::string materialName = "";
        strncpy(m_primInteractMaterialName, materialName.c_str(),80);
        m_primInteractFlag = -1;
        m_primInteractX0 = -999.;
        m_primInteractY0 = -999.;
        m_primInteractZ0 = -999.;
        m_primInteractEleXDir = -999.;
        m_primInteractEleYDir = -999.;
        m_primInteractEleZDir = -999.;
        m_primInteractEleEne  = -999.;
        m_primInteractPosXDir = -999.;
        m_primInteractPosYDir = -999.;
        m_primInteractPosZDir = -999.;
        m_primInteractPosEne  = -999.;
        double MaxDocaDght = 0.0001;
        if (m_primType == 22){ // we start with a photon 
          
          HepPoint3D primaryPosition              = primary->initialPosition() ;
          CLHEP::HepLorentzVector primaryMomentum = primary->initialFourMomentum() ;

          // Need pointers to the Ele/Pos daughters
          const Event::McParticle* daughterEle = 0;
          const Event::McParticle* daughterPos = 0;
          double daughterElePrimaryDoca = 100000.;
          double daughterPosPrimaryDoca = 100000.;

          // Need the Ray(Point, Vector) of the primary particle for the doca
          Point  x1 = Point(primaryPosition.x(), primaryPosition.y(), primaryPosition.z() );
          Vector t1 = Vector(primaryMomentum.x(), primaryMomentum.y(), primaryMomentum.z()).unit();
          Doca primaryRay(x1, t1);

          // Loop through the daughters 
          const SmartRefVector<Event::McParticle>& daughterList = primary->daughterList();

          for(SmartRefVector<Event::McParticle>::const_iterator daughterIter = daughterList.begin();
              daughterIter != daughterList.end(); daughterIter++)
            {
              const Event::McParticle* daughter = *daughterIter;

              // get the electron with initial position closest to the initial photon direction 
              if (daughter->particleProperty() == 11){
                Point  pp(daughter->initialPosition().x(), daughter->initialPosition().y(), daughter->initialPosition().z() );
                double doca = primaryRay.docaOfPoint( pp );
                if (doca < daughterElePrimaryDoca)
                  {
                    daughterEle            = daughter;
                    daughterElePrimaryDoca = doca;
                  }
              }
          
              // get the positron with initial position closest to the initial photon direction 
              if (daughter->particleProperty() == -11){
                Point  pp(daughter->initialPosition().x(), daughter->initialPosition().y(), daughter->initialPosition().z() );
                double doca = primaryRay.docaOfPoint( pp );
                if (doca < daughterPosPrimaryDoca)
                  {
                    daughterPos            = daughter;
                    daughterPosPrimaryDoca = doca;
                  }
              }
            }// end daughter loop
            

          // if they both are from pair prod [and start form same position] and (DOCA < MaxDocaDght): Pair production found!
          // if electron is from compton and (DOCA < MaxDocaDght) [and no positron] : Compton found!
          // else: we don't care! 
          std::string convString  = "conv" ;
          std::string comptString = "compt";
          std::string processEle = daughterEle ? daughterEle->getProcess() : "";
          std::string processPos = daughterPos ? daughterPos->getProcess() : "";

          if ( (processEle==convString) && (processPos==convString) && daughterElePrimaryDoca< MaxDocaDght && daughterPosPrimaryDoca< MaxDocaDght){
            // if pair production use photon final position
            m_primInteractFlag = 0;
            m_primInteractX0 = primary->finalPosition().x();
            m_primInteractY0 = primary->finalPosition().y();
            m_primInteractZ0 = primary->finalPosition().z();
            // set primary electron momentum
            m_primInteractEleXDir = daughterEle->initialFourMomentum().vect().unit().x();
            m_primInteractEleYDir = daughterEle->initialFourMomentum().vect().unit().y();
            m_primInteractEleZDir = daughterEle->initialFourMomentum().vect().unit().z();
            m_primInteractEleEne  = daughterEle->initialFourMomentum().t();
            // set primary positron momentum
            m_primInteractPosXDir = daughterPos->initialFourMomentum().vect().unit().x();
            m_primInteractPosYDir = daughterPos->initialFourMomentum().vect().unit().y();
            m_primInteractPosZDir = daughterPos->initialFourMomentum().vect().unit().z();
            m_primInteractPosEne  = daughterPos->initialFourMomentum().t();

          }
          else if ((processEle==comptString) && daughterElePrimaryDoca< MaxDocaDght ){
            // if compton, use electron initial position
            m_primInteractFlag = 1;
            m_primInteractX0 = daughterEle->initialPosition().x();
            m_primInteractY0 = daughterEle->initialPosition().y();
            m_primInteractZ0 = daughterEle->initialPosition().z();
            // set primary electron momentum
            m_primInteractEleXDir = daughterEle->initialFourMomentum().vect().unit().x();
            m_primInteractEleYDir = daughterEle->initialFourMomentum().vect().unit().y();
            m_primInteractEleZDir = daughterEle->initialFourMomentum().vect().unit().z();
            m_primInteractEleEne  = daughterEle->initialFourMomentum().t();
          }
          else {
            // if other cases production use photon final position
            m_primInteractFlag = 2;
            m_primInteractX0 = primary->finalPosition().x();
            m_primInteractY0 = primary->finalPosition().y();
            m_primInteractZ0 = primary->finalPosition().z();
          }

          // now use the propagation tool to access material Id and geometry. (thanks to tracy - jira LPATE-81)
          Point  pp(m_primInteractX0, m_primInteractY0, m_primInteractZ0);
          m_G4PropTool->setStepStart(pp, t1);
          m_primInteractMaterialId = m_G4PropTool->getMaterialIndex(0);
          materialName = m_G4PropTool->getMaterialName(0);
          strncpy(m_primInteractMaterialName, materialName.c_str(),80);
          if (m_primInteractMaterialName == "") memset(m_primInteractMaterialName, '_', 80);

        } // end if photon.
    }

    // Check to see that **something** happened
    m_primInteracted = GetInteracted();

    // Clean up after ourselves
    if (m_partToTrajTab)    {delete m_partToTrajTab;    m_partToTrajTab    = 0;}
    if (m_pointToPosHitTab) {delete m_pointToPosHitTab; m_pointToPosHitTab = 0;}
    if (m_pointToIntHitTab) {delete m_pointToIntHitTab; m_pointToIntHitTab = 0;}
    if (m_clusToPosHitTab)  {delete m_clusToPosHitTab;  m_clusToPosHitTab  = 0;}
    
    log << MSG::DEBUG << " returning. " << endreq;

    return sc;
}


void McTkrHitValsTool::CntTotalHits(const Event::McParticle* /* primary */)
{
    m_totalTkrHits     = 0;
    m_totalPrimTkrHits = 0;

    // If there is no cluster table then there are no hits (by definition)
    if (m_clusToPosHitTab && !m_clusToPosHitTab->getAllRelations()->empty())
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
            m_totalTkrHits += idToHitIter->second.size();
        }

        // Go through and count total number of hits from the primary and its descendants
        for (std::set<int>::iterator trkIdIter = primTrkIdSet.begin(); trkIdIter != primTrkIdSet.end(); trkIdIter++)
        {
            m_totalPrimTkrHits += idToHitMap[*trkIdIter].size();
        }

        // And now the total number of clusters from the primary and its descendants
        m_totalPrimClusHits = primClusSet.size();
    }

    return;
}

int McTkrHitValsTool::GetAcdHits(const Event::McParticle* particle)
{
    int nAcdHits = 0;

    // If there is no point to position hit table then there are no hits
    if (m_pointToPosHitTab)
    {
        // Look up trajectory for this particle
        Event::McPartToTrajectoryVec trajVec = m_partToTrajTab->getRelByFirst(particle);
        Event::McTrajectory* mcTraj          = (*(trajVec.begin()))->getSecond();

        // Set up to loop through Trajectory points
        std::vector<Event::McTrajectoryPoint*> points = mcTraj->getPoints();
        std::vector<Event::McTrajectoryPoint*>::const_iterator pointIter;

        // Keep track of number of McPositionHits in each ACD object
        std::map<idents::AcdId, std::vector<const Event::McPositionHit*> > acdToHitMap;
        acdToHitMap.clear();

        // Loop through all the McTrajectoryPoints on this track
        for(pointIter = points.begin(); pointIter != points.end(); pointIter++) 
        {
            // De-reference the McTrajectoryPoint pointer
            Event::McTrajectoryPoint* mcPoint = *pointIter;

            // Looking for ACD hits only
            if (!mcPoint->getVolumeID().isAcd()) continue;

            // Of course, that only gives us the general volume... must do rest inside a try-catch
            try
            {
                idents::AcdId acdId(mcPoint->getVolumeID());

                // Look for McTrajectoryPoint <--> McPositionHit
                Event::McPointToPosHitVec posRelVec = m_pointToPosHitTab->getRelByFirst(mcPoint);

                // Did this trajectory point leave a mark?
                if (!posRelVec.empty())
                {
                    // Can only be one McPositionHit per McTrajectoryPoint... by definition
                    Event::McPointToPosHitRel* hitRel = *(posRelVec.begin());

                    Event::McPositionHit* mcPosHit = hitRel->getSecond();

                    // Add to local id to posHit map for counting
                    acdToHitMap[acdId].push_back(mcPosHit);

                    // Update global id to energy map
                    // For ACD, McPositionHits are unique for each McParticle so not double counting here
                    m_acdIdToEnergyMap[acdId] += mcPosHit->depositedEnergy();
                }
            }
            catch(...)
            {
                continue;
            }
        }

        // ok, the size of trackClusters is, then, the number of hits this track left
        nAcdHits = acdToHitMap.size();
    }

    return nAcdHits;
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

            // Tracker only
            if (!mcPoint->getVolumeID().isTkr()) continue;

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

int McTkrHitValsTool::GetCalHits(const Event::McParticle* particle)
{
    int nCalHits = 0;

    // If there is no cluster table then there are no hits
    if (m_pointToIntHitTab)
    {
        // Look up trajectory for this particle
        Event::McPartToTrajectoryVec trajVec = m_partToTrajTab->getRelByFirst(particle);
        Event::McTrajectory* mcTraj          = (*(trajVec.begin()))->getSecond();

        // Set up to loop through Trajectory points
        std::vector<Event::McTrajectoryPoint*> points = mcTraj->getPoints();
        std::vector<Event::McTrajectoryPoint*>::const_iterator pointIter;

        // Keep track of number of McIntegratingHits in each Xtal object
        std::map<idents::CalXtalId, std::vector<const Event::McIntegratingHit*> > xtalToIntMap;
        xtalToIntMap.clear();

        // Loop through all the McTrajectoryPoints on this track
        for(pointIter = points.begin(); pointIter != points.end(); pointIter++) 
        {
            // De-reference the McTrajectoryPoint pointer
            Event::McTrajectoryPoint* mcPoint = *pointIter;

            // Looking for Cal hits only
            if (!mcPoint->getVolumeID().isCal()) continue;

            // And, like with ACD, this only gets us close, rest in a try-catch block
            try
            {
                idents::CalXtalId calXtalId(mcPoint->getVolumeID());

                // Look for McTrajectoryPoint <--> McIntegratingHit
                Event::McPointToIntHitVec intRelVec = m_pointToIntHitTab->getRelByFirst(mcPoint);

                // Did this trajectory point leave a mark?
                if (!intRelVec.empty())
                {
                    // Can only be one McIntegratingHit per McTrajectoryPoint... by definition
                    Event::McPointToIntHitRel* hitRel = *(intRelVec.begin());

                    Event::McIntegratingHit* mcIntHit = hitRel->getSecond();

                    xtalToIntMap[calXtalId].push_back(mcIntHit);
                }
            }
            catch(...)
            {
                continue;
            }
        }

        // For Calorimeter, McIntegratingHits are unique to xtalIds, must loop here and update
        for(std::map<idents::CalXtalId, std::vector<const Event::McIntegratingHit*> >::iterator mapIter  = xtalToIntMap.begin();
                                                                                                mapIter != xtalToIntMap.end();
                                                                                                mapIter++)
        {
            if (m_xtalIdToEnergyMap.find(mapIter->first) == m_xtalIdToEnergyMap.end())
            {
                m_xtalIdToEnergyMap[mapIter->first] = mapIter->second.front()->totalEnergy();
            }
        }

        // Set number of hit crystals by simulated event
        nCalHits = xtalToIntMap.size();
    }

    return nCalHits;
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
            if (!volId.isTkr()) continue;

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

int  McTkrHitValsTool::GetInteracted()
{
    int interacted = 0;

    // First check McPositionHits to see if any were created by the simulated
    // and are not from the Overlay
    SmartDataPtr<Event::McPositionHitCol> posHitCol(m_pEventSvc, EventModel::MC::McPositionHitCol);

    if (posHitCol)
    {
        // Loop through the McPositionHit collection
        for (Event::McPositionHitCol::iterator posHitIter = posHitCol->begin(); posHitIter != posHitCol->end(); posHitIter++)
        {
            Event::McPositionHit* posHit = *posHitIter;

            // If the overlay bit is NOT set then this is a hit from simulation
            if ((posHit->getPackedFlags() & (Event::McPositionHit::overlayHit | Event::McPositionHit::simulation))
                                          != Event::McPositionHit::overlayHit)
            {
                interacted = 1;
                break;
            }
        }
    }

    if (!interacted)
    {
        // Now look at the Calorimeter
        SmartDataPtr<Event::McIntegratingHitCol> intHitCol(m_pEventSvc, EventModel::MC::McIntegratingHitCol);

        if (intHitCol)
        {
            // Loop through the McIntegratingHit collection
            for (Event::McIntegratingHitCol::iterator intHitIter = intHitCol->begin(); intHitIter != intHitCol->end(); intHitIter++)
            {
                Event::McIntegratingHit* intHit = *intHitIter;

                // Ditto as above, if overlay bit not set this is hit from simulation
                if ((intHit->getPackedFlags() & (Event::McIntegratingHit::overlayHit | Event::McIntegratingHit::simulation))
                                              != Event::McIntegratingHit::overlayHit)
                {
                    interacted = 1;
                    break;
                }
            }
        }
    }

    return interacted;
}


void McTkrHitValsTool::countTruncatedPlanes(Event::McParticle* primary)
{
    // Apparently no way around the issue of having to set the array of active distances to guard values...
    for(int idx = 0; idx < 36; idx++) m_truncActDist[idx] = -m_towerPitch;

    // Set up to look for truncated planes
    SmartDataPtr<Event::McPartToTrajectoryTabList> mcPartToTraj(m_pEventSvc, "/Event/MC/McPartToTrajectory");

    if (mcPartToTraj)
    {
        // Recover the full relational table
        Event::McPartToTrajectoryTab mcPartToTrajTab(mcPartToTraj.ptr());
    
        // Recover the McPositionHits to trajectory points relational tables
        SmartDataPtr<Event::McPointToPosHitTabList> mcPointToPosHitList(m_pEventSvc, "/Event/MC/McPointToPosHit");

        if (!mcPointToPosHitList)
        { 
            mcPointToPosHitList = new Event::RelationList<Event::McTrajectoryPoint, Event::McPositionHit>;
        }
        Event::McPointToPosHitTab mcPosHitTab(mcPointToPosHitList.ptr());

        // Get the truncation information
        SmartDataPtr<Event::TkrTruncationInfo> truncInfo( m_pEventSvc, EventModel::TkrRecon::TkrTruncationInfo );

        // Ultimate goal is map between planes and active distance to truncated regions
        std::map<int, bool>    planeIsTruncatedMap;
        std::map<int, double>  planeToActDistMap;

        // Get the daughter list we will loop over 
        SmartRefVector<Event::McParticle> daughters = primary->daughterList();

        // Loop over daughters of the primary. In theory, for gammas, there should be two which can be
        // making tracks. We'll look at each independently.
        for(SmartRefVector<Event::McParticle>::iterator mcPartItr = daughters.begin(); mcPartItr != daughters.end(); mcPartItr++)
        {
            Event::McParticle*           mcPart  = *mcPartItr;
            Event::McPartToTrajectoryVec trajVec = mcPartToTrajTab.getRelByFirst(mcPart);

            // Goal is to build a map between planes crossed and TkrTruncatedPlane info
            std::map<Event::SortId, Event::TkrTruncatedPlane*> planeToTruncInfoMap;
            std::map<Event::SortId, CLHEP::Hep3Vector>         planeToPositionMap;

            // This first task is to go through the trajectory associated with this daughter particle
            // and store away its position in the sensitive volumes and keep a pointer to the truncation
            // information for that plane
            if (!trajVec.empty())
            {
                Event::McTrajectory* traj = trajVec.front()->getSecond();

                std::vector<Event::McTrajectoryPoint*> points = traj->getPoints();
                std::vector<Event::McTrajectoryPoint*>::const_iterator pit;

                // Loop the trajectory points to get positions in the sensitive planes
                for(pit = points.begin(); pit != points.end(); pit++) 
                {
                    Event::McTrajectoryPoint* point = *pit;

                    // If in the tracker, next is to make sure we are in a sensitive volume
                    // Test for this is that an McPositionHit exists (is there a better way?)
                    Event::McPointToPosHitVec relVec = mcPosHitTab.getRelByFirst(point);

                    if (relVec.empty()) continue;

                    idents::VolumeIdentifier  volId = relVec.front()->getSecond()->volumeID();

                    // If not in the tracker then nothing to do here
                    if (!volId.isTkr()) continue;

                    // Convert the volume identifier into 
                    idents::TkrId tkrId(volId);

                    // Can this happen?
                    if (!tkrId.hasTray() || !tkrId.hasBotTop() || !tkrId.hasView())
                    {
                        continue;
                    }

                    // Convert this into a SortId
                    Event::SortId sortId(idents::TowerId(tkrId.getTowerX(), tkrId.getTowerY()).id(),
                                         tkrId.getTray(), tkrId.getBotTop(), tkrId.getView());

                    // We can have multiple points within a given sensitive plane as the particle
                    // interacts as it enters and then exists the silicon. For now we are only
                    // interested in the position as it enters the silicon - presumed from above
                    std::map<Event::SortId, CLHEP::Hep3Vector>::iterator planeToPosItr = planeToPositionMap.find(sortId);

                    if (planeToPosItr != planeToPositionMap.end() && planeToPosItr->second.z() > point->getPoint().z())
                    {
                        continue;
                    }

                    // Null pointer for default situation
                    Event::TkrTruncatedPlane* truncPlane = 0;

                    // If truncation information then see what we can find
                    if (truncInfo)
                    {
                        Event::TkrTruncationInfo::TkrTruncationMap::iterator planeItr = truncInfo->getTruncationMap()->find(sortId);

                        if (planeItr != truncInfo->getTruncationMap()->end()) truncPlane = &planeItr->second;
                    }

                    // Store entry in map
                    planeToTruncInfoMap[sortId] = truncPlane;
                    planeToPositionMap[sortId]  = point->getPoint();
                }
            }

            // Once done looping over daughters then try to decode the information we have
            // The plan is to run through the map and reorganize by planes, but this time
            // only keeping the "active distance" information
            for(std::map<Event::SortId, CLHEP::Hep3Vector>::iterator idItr  = planeToPositionMap.begin(); 
                                                                     idItr != planeToPositionMap.end(); 
                                                                     idItr++)
            {
                Event::SortId     sortId   = idItr->first;
                CLHEP::Hep3Vector position = idItr->second;

                // Convert back to tracker id to get plane
                idents::TowerId tower(sortId.getTower());
                idents::TkrId   tkrId(tower.ix(), tower.iy(), sortId.getTray(), sortId.getFace(), sortId.getView());

                int plane = m_tkrGeom->getPlane(tkrId);

                // Make record for all time
                std::map<int, bool>::iterator planeIsTruncatedItr = planeIsTruncatedMap.find(plane);

                if (planeIsTruncatedItr == planeIsTruncatedMap.end())
                {
                    planeIsTruncatedMap[plane] = false;
                    planeToActDistMap[plane]   = -10000.;
                }

                // Recover the pointer to the truncated plane information
                Event::TkrTruncatedPlane* truncatedPlane = planeToTruncInfoMap[sortId];

                if (truncatedPlane)
                {
                    const int status = truncatedPlane->getStatus();

                    if (!status) continue;

                    const Event::intVector&   numStrips  = truncatedPlane->getStripCount();
                    const Event::floatVector& localX     = truncatedPlane->getLocalX();
                    double                    z          = truncatedPlane->getPlaneZ();
                    double                    truncScale = 4;
                    double                    dz         = 0.5 * m_siThickness * truncScale;
                    
                    idents::TowerId towerId(sortId.getTower());
                    
                    int    xTower  = towerId.ix();
                    int    yTower  = towerId.iy();
                    double offSetX = m_towerPitch * (xTower - 0.5 * (m_numXTowers - 1));
                    double offSetY = m_towerPitch * (yTower - 0.5 * (m_numYTowers - 1));

                    double xLow  = localX[3];
                    double xHigh = localX[3];

                    if((status & Event::TkrTruncatedPlane::END0SET)>0)                   xLow  = localX[0];
                    if((status & Event::TkrTruncatedPlane::RC1SET)>0)                    xHigh = localX[1];
                    if((status & Event::TkrTruncatedPlane::CC1SET)>0 && numStrips[1]==0) xHigh = localX[2];

                    double activeDist = planeToActDistMap[plane];

                    double xPointPos = position.x() - offSetX;
                    double yPointPos = position.y() - offSetY;

                    if (sortId.getView() == idents::TkrId::eMeasureY)
                        std::swap(xPointPos, yPointPos);

                    if(xHigh > xLow) 
                    {
                        double x  = 0.5 * (xLow + xHigh);
                        double dx = 0.5 * (xHigh - xLow);
                        double y  = 0.0;
                        double dy = 0.5 * m_stripLength;

                        double dMeasPos   = xPointPos - x;
                        double dUnMeasPos = yPointPos - y;

                        double activeDistMeas   = dx - fabs(dMeasPos);
                        double activeDistUnMeas = dy - fabs(dUnMeasPos);

                        // Now, want a specific value here... 
                        // Assume we are IN a truncated region, want to return the distance in
                        // the measured direction. Also, if both our "out" then also return the
                        // same distance.
                        activeDist = activeDistMeas;

                        // If measured is in but unmeasured is out, return that one...
                        if (activeDistMeas >=0. && activeDistUnMeas < 0.) 
                        {
                            activeDist = activeDistUnMeas;
                        }
                    }
                    if((status&Event::TkrTruncatedPlane::CC1SET)>0 && numStrips[1]>0){
                        xLow  = localX[2];
                        xHigh = 0.5 * m_activeWidth;

                        double x  = 0.5 * (xLow + xHigh);
                        double dx = 0.5 * (xHigh - xLow);
                        double y  = 0.0;
                        double dy = 0.5 * m_stripLength;

                        double dMeasPos   = xPointPos - x;
                        double dUnMeasPos = yPointPos - y;

                        double activeDistMeas   = dx - fabs(dMeasPos);
                        double activeDistUnMeas = dy - fabs(dUnMeasPos);

                        // Now, want a specific value here... 
                        // Assume we are IN a truncated region, want to return the distance in
                        // the measured direction. Also, if both our "out" then also return the
                        // same distance.
                        double localActiveDist = activeDistMeas;

                        // If measured is in but unmeasured is out, return that one...
                        if (activeDistMeas >=0. && activeDistUnMeas < 0.) 
                        {
                            localActiveDist = activeDistUnMeas;
                        }

                        // Combine this result with the possible previous result (though we believe this can't happen)
                        activeDist = std::max(activeDist, localActiveDist);
                    }

                    // Store away
                    if (activeDist > planeToActDistMap[plane])
                    {
                        planeIsTruncatedMap[plane] = true;           // The plane is truncated
                        planeToActDistMap[plane]   = activeDist;     // Active distance to the hit
                    }
                }
            }
        }

        // Ok, now extract some useful information, but only if there is actually something there!
        if (!planeIsTruncatedMap.empty())
        {
            // First is easy since map arranges bottom plane first in list
            m_bottomPlane = planeIsTruncatedMap.begin()->first;

            // The number of planes crossed is the number of map entries 
            // which is NOT the same as the total number of planes from top to bottom!
            m_planesCrossed = planeIsTruncatedMap.size();

            // Loop through entries in the map, this will be in ascending plane order
            for(std::map<int, bool>::iterator planeTruncItr = planeIsTruncatedMap.begin(); planeTruncItr != planeIsTruncatedMap.end(); planeTruncItr++)
            {
                int plane = planeTruncItr->first;

                if (m_topPlane < plane) m_topPlane = plane;

                if (planeTruncItr->second)
                {
                    m_nTruncatedPlanes++;

                    m_truncActDist[plane] = planeToActDistMap[plane];

                    if (planeToActDistMap[plane] >= 0) m_nTruncatedLost++;
                }
            }

            int success = 10000;
        }
    }

    return;
}
