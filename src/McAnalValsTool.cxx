/** @file McAnalValsTool.cxx
    @brief declartion, implementaion of the class UserAlg

    $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/McAnalValsTool.cxx,v 1.6 2004/02/20 02:30:01 lsrea Exp $
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
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "GlastSvc/MonteCarlo/IMcGetEventInfoTool.h"
#include "GlastSvc/MonteCarlo/IMcGetTrackInfoTool.h"

#include <algorithm>

// to activate debug vars
// #define DEBUG

/*! @class McAnalValsTool
@brief calculates Monte Carlo values

@authors Bill Atwood, Leon Rochester
*/

class McAnalValsTool : public ValBase
{
public:
    
    McAnalValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~McAnalValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    /// For determing MC tuple values
    void       calcMcAngleInfo(const Event::McParticle* mcPart, const Event::McParticle* mcGamma,
                               double& ang2Gam, double& cls2Gam, double& mcTrkRms, double& prtType,
                               double& prtEnergy, double& prtCosDirX, double& prtCosDirY, double& prtCosDirZ,
                               double& prtNHits, double& prtNClstrs, double& prtNGaps, double& prt1stGapSz,
                               double& prtNHits2Gp);

    void       calcMcEnergyInfo(const Event::McParticle* mcPart,
                                double& lastHitE, double& radLossE, double& dltaRayE, double& nBrems,
                                double& nDeltas, double&nDeltaHt, double& aveRange, double& maxRange);

    //Pure MC Tuple Items
    double       m_numCalls;

    // Variables for the Primary track (the gamma in gamma events)
    double       m_prmEnergy;         // Energy of primary at creation
    double       m_prmPosX;           // Start position X
    double       m_prmPosY;           // Start position Y
    double       m_prmPosZ;           // Start position Z
    double       m_prmCosDirX;        // Start direction cosine X
    double       m_prmCosDirY;        // Start direction cosine Y
    double       m_prmCosDirZ;        // Start direction cosine Z
    double       m_prmDecEne;         // Energy of primary at decay/stop
    double       m_prmDecPosX;        // Decay position X
    double       m_prmDecPosY;        // Decay position Y
    double       m_prmDecPosZ;        // Decay position Z
    double       m_prmDecCosX;        // Decay direction cosine X
    double       m_prmDecCosY;        // Decay direction cosine Y
    double       m_prmDecCosZ;        // Decay direction cosine Z
    double       m_prmNDghtrs;        // Number of daughters in decay
    double       m_prmDecCode;        // Decay classification (see McEventStructure.h)

    double       m_prmDecLayer;       // Tray number for decay
    double       m_prmDecInCnv;       // If decay in tracker tray, material ID

    double       m_prmNumSecndry;     // Number of "secondaries" (McTracks)
    double       m_prmNumAsscted;     // Number of "associated" particles
    double       m_prmTrkPattern;     // Monte Carlo track hit pattern
    double       m_prmMcAngle;        // Angle between primary and resultant from secondaries
    double       m_prmClsAngle;       // Same, but using cluster positions

    // Variables for the "first" secondary in decay (or further info for primary if charged event)
    double       m_scd1Type;          // Type of particle
    double       m_scd1FirstLyr;      // First hit layer of track
    double       m_scd1LastLyr;       // Last hit layer of track
    double       m_scd1Energy;        // Energy of particle at creation
    double       m_scd1CosDirX;       // Direction cosine X
    double       m_scd1CosDirY;       // Direction cosine Y
    double       m_scd1CosDirZ;       // Direction cosine Z
    double       m_scd1NHits;         // Number of McSiLayerHits made by particle
    double       m_scd1NClstrs;       // Number of clusters associated to particle
    double       m_scd1NGaps;         // Number of gaps in layers particle makes in tracker
    double       m_scd11stGapSz;      // Size of the first gap
    double       m_scd1NHits2Gp;      // Number of hits before first gap
    double       m_scd1McTrkRms;      // RMS of the track
    double       m_scd1Ang2Gam;       // Angle the particle makes to the primary
    double       m_scd1Cls2Gam;       // Same but using clusters
    double       m_scd1ELastHit;      // Energy of the track at the last hit in the tracker
    double       m_scd1RadELoss;      // Energy radiated by the particle from birth to last hit
    double       m_scd1DeltaRay;      // Energy lost to delta rays from birth to last hit
    double       m_scd1NBrems;        // Number of bremstrahlung photons produced from birth to last hit
    double       m_scd1NDeltas;       // Number of delta rays produced from birth to last hit
    double       m_scd1NDeltaHt;      // Number of above which create McPositionHits in the tracker
    double       m_scd1AveRange;      // Average range of all deltas produced from birth to last hit
    double       m_scd1MaxRange;      // Maximum range of deltas produced from birth to last hit

    // Variables for the "second" secondary - see above definitions
    double       m_scd2Type;
    double       m_scd2FirstLyr;      
    double       m_scd2LastLyr;       
    double       m_scd2Energy;
    double       m_scd2CosDirX;
    double       m_scd2CosDirY;
    double       m_scd2CosDirZ;
    double       m_scd2NHits;
    double       m_scd2NClstrs;
    double       m_scd2NGaps;
    double       m_scd21stGapSz;
    double       m_scd2NHits2Gp;
    double       m_scd2McTrkRms;
    double       m_scd2Ang2Gam;
    double       m_scd2Cls2Gam;
    double       m_scd2ELastHit;
    double       m_scd2RadELoss;
    double       m_scd2DeltaRay;
    double       m_scd2NBrems;
    double       m_scd2NDeltas;
    double       m_scd2NDeltaHt;
    double       m_scd2AveRange;
    double       m_scd2MaxRange;

    // to decode the particle charge
    IParticlePropertySvc*  m_ppsvc;

    // Keep track of tool for mc tracks
    IMcGetEventInfoTool*   m_mcEvent;
    IMcGetTrackInfoTool*   m_mcTracks;
};

// Static factory for instantiation of algtool objects
static ToolFactory<McAnalValsTool> s_factory;
const IToolFactory& McAnalValsToolFactory = s_factory;

// Standard Constructor
McAnalValsTool::McAnalValsTool(const std::string& type, 
                               const std::string& name, 
                               const IInterface* parent)
                               : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode McAnalValsTool::initialize()
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

    // TO DO here: gracefully return if tools not located, set up to NOT run the tool
    m_mcEvent = 0;
    sc = toolSvc()->retrieveTool("McGetEventInfoTool", m_mcEvent);
    if (sc.isFailure()) {
        log << MSG::INFO << " McGetEventInfoTool not found" << endreq;
        log << MSG::INFO << " Will not generate McAnalVals" << endreq;
        return StatusCode::SUCCESS;
    }

    m_mcTracks = 0;
    sc = toolSvc()->retrieveTool("McGetTrackInfoTool", m_mcTracks);
    if (sc.isFailure()) {
        log << MSG::INFO << " McGetTrackInfoTool not found!" << endreq;
        log << MSG::INFO << " Will not generate McAnalVals" << endreq;
        return StatusCode::SUCCESS;
    }
    
    // load up the map

	addItem("McaNumCalls",       &m_numCalls);
	addItem("McaPrmEnergy",      &m_prmEnergy);
	addItem("McaPrmDecEne",      &m_prmDecEne);
	addItem("McaDecPosX",        &m_prmPosX);
	addItem("McaDecPosY",        &m_prmPosY);
	addItem("McaDecPosZ",        &m_prmPosZ);
	addItem("McaPrmDecPosX",     &m_prmDecPosX);
	addItem("McaPrmDecPosY",     &m_prmDecPosY);
	addItem("McaPrmDecPosZ",     &m_prmDecPosZ);
	addItem("McaPrmCosDirX",     &m_prmCosDirX);
	addItem("McaPrmCosDirY",     &m_prmCosDirY);
	addItem("McaPrmCosDirZ",     &m_prmCosDirZ);
	addItem("McaPrmDecCosX",     &m_prmDecCosX);
	addItem("McaPrmDecCosY",     &m_prmDecCosY);
	addItem("McaPrmDecCosZ",     &m_prmDecCosZ);
	addItem("McaPrmNDghtrs",     &m_prmNDghtrs);
	addItem("McaPrmDecCode",     &m_prmDecCode);


    addItem("McaPrmNumSecndry",  &m_prmNumSecndry);
    addItem("McaPrmNumAsscted",  &m_prmNumAsscted);
    addItem("McaPrmTrkPattern",  &m_prmTrkPattern);
    addItem("McaPrmMcAngle",     &m_prmMcAngle);
    addItem("McaPrmClsAngle",    &m_prmClsAngle);

    addItem("McaPrmDecLayer",    &m_prmDecLayer);
    addItem("McaPrmDecInCnv",    &m_prmDecInCnv);

    addItem("McaScd1PartType",   &m_scd1Type);
    addItem("McaScd1FirstLyr",   &m_scd1FirstLyr);
    addItem("McaScd1LastLyr",    &m_scd1LastLyr);
    addItem("McaScd1Energy",     &m_scd1Energy);
    addItem("McaScd1CosDirX",    &m_scd1CosDirX);
    addItem("McaScd1CosDirY",    &m_scd1CosDirY);
    addItem("McaScd1CosDirZ",    &m_scd1CosDirZ);
    addItem("McaScd1NHits",      &m_scd1NHits);
    addItem("McaScd1NClstrs",    &m_scd1NClstrs);
    addItem("McaScd1NGaps",      &m_scd1NGaps);
    addItem("McaScd11stGapSz",   &m_scd11stGapSz);
    addItem("McaScd1NHits2Gp",   &m_scd1NHits2Gp);
    addItem("McaScd1McTrkRms",   &m_scd1McTrkRms);
    addItem("McaScd1Ang2Gam",    &m_scd1Ang2Gam);
    addItem("McaScd1Cls2Gam",    &m_scd1Cls2Gam);
    addItem("McaScd1ELastHit",   &m_scd1ELastHit);
    addItem("McaScd1RadEloss",   &m_scd1RadELoss);
    addItem("McaScd1DeltaRay",   &m_scd1DeltaRay);
    addItem("McaScd1NBrems",     &m_scd1NBrems);
    addItem("McaScd1NDeltas",    &m_scd1NDeltas);
    addItem("McaScd1NDeltaHt",   &m_scd1NDeltaHt);
    addItem("McaScd1AveRange",   &m_scd1AveRange);
    addItem("McaScd1MaxRange",   &m_scd1MaxRange);

    addItem("McaScd2PartType",   &m_scd2Type);
    addItem("McaScd2FirstLyr",   &m_scd2FirstLyr);
    addItem("McaScd2LastLyr",    &m_scd2LastLyr);
    addItem("McaScd2Energy",     &m_scd2Energy);
    addItem("McaScd2CosDirX",    &m_scd2CosDirX);
    addItem("McaScd2CosDirY",    &m_scd2CosDirY);
    addItem("McaScd2CosDirZ",    &m_scd2CosDirZ);
    addItem("McaScd2NHits",      &m_scd2NHits);
    addItem("McaScd2NClstrs",    &m_scd2NClstrs);
    addItem("McaScd2NGaps",      &m_scd2NGaps);
    addItem("McaScd21stGapSz",   &m_scd21stGapSz);
    addItem("McaScd2NHits2Gp",   &m_scd2NHits2Gp);
    addItem("McaScd2McTrkRms",   &m_scd2McTrkRms);
    addItem("McaScd2Ang2Gam",    &m_scd2Ang2Gam);
    addItem("McaScd2Cls2Gam",    &m_scd2Cls2Gam);
    addItem("McaScd2ELastHit",   &m_scd2ELastHit);
    addItem("McaScd2RadEloss",   &m_scd2RadELoss);
    addItem("McaScd2DeltaRay",   &m_scd2DeltaRay);
    addItem("McaScd2NBrems",     &m_scd2NBrems);
    addItem("McaScd2NDeltas",    &m_scd2NDeltas);
    addItem("McaScd2NDeltaHt",   &m_scd2NDeltaHt);
    addItem("McaScd2AveRange",   &m_scd2AveRange);
    addItem("McaScd2MaxRange",   &m_scd2MaxRange);

    zeroVals();
    
    return sc;
}


StatusCode McAnalValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream  log( msgSvc(), name() );

    if (m_mcEvent==0 || m_mcTracks==0) return sc;

    // Retrieving pointers from the TDS 
    SmartDataPtr<Event::EventHeader>   header(m_pEventSvc,    EventModel::EventHeader);
    SmartDataPtr<Event::MCEvent>       mcheader(m_pEventSvc,  EventModel::MC::Event);
    SmartDataPtr<Event::McParticleCol> particles(m_pEventSvc, EventModel::MC::McParticleCol);

    double t = header->time();
    log << MSG::DEBUG << "Event time: " << t << endreq;;

    if (m_mcEvent)
    {
        //int numTracksTotal = m_mcEvent->getNumMcTracks();
        int classifyBits   = m_mcEvent->getClassificationBits();

        // Pointers to the primary particle and (eventually) its secondaries
        const Event::McParticle* mcPart   = m_mcEvent->getPrimaryParticle();
        const Event::McParticle* mcMain   = 0;
        const Event::McParticle* mcSecond = 0;

        idents::VolumeIdentifier curVolId    = const_cast<Event::McParticle*>(mcPart)->getFinalId();
        HepLorentzVector         primary4mom = mcPart->initialFourMomentum();

        // Retrieve information about the primary particle
        m_prmEnergy     = primary4mom.e();
        m_prmDecEne     = mcPart->finalFourMomentum().e();
        m_prmPosX       = mcPart->initialPosition().x();
        m_prmPosY       = mcPart->initialPosition().y();
        m_prmPosZ       = mcPart->initialPosition().z();
        m_prmCosDirX    = primary4mom.vect().unit().x();
        m_prmCosDirY    = primary4mom.vect().unit().y();
        m_prmCosDirZ    = primary4mom.vect().unit().z();
        m_prmDecPosX    = mcPart->finalPosition().x();
        m_prmDecPosY    = mcPart->finalPosition().y();
        m_prmDecPosZ    = mcPart->finalPosition().z();
        m_prmDecCosX    = mcPart->finalFourMomentum().vect().unit().x();
        m_prmDecCosY    = mcPart->finalFourMomentum().vect().unit().y();
        m_prmDecCosZ    = mcPart->finalFourMomentum().vect().unit().z();
        m_prmDecCode    = classifyBits;
        m_prmNDghtrs    = mcPart->daughterList().size();
        m_prmNumSecndry = m_mcEvent->getNumSecondaries();
        m_prmNumAsscted = m_mcEvent->getNumAssociated();

        // If charged particle then get base information from this path...
        if (classifyBits & Event::McEventStructure::CHARGED)
        {
            // Get the initial volume id
            curVolId = const_cast<Event::McParticle*>(mcPart)->getInitialId();

            // For charged tracks, the main secondary is the track itself
            mcMain   = mcPart;
        }
        // If gamma then this path
        else if (classifyBits & Event::McEventStructure::GAMMA)
        {
            primary4mom = mcPart->finalFourMomentum();

            // This should mean that the gamma converted in the tracker (or, unfortunately, the ACD)
            if (m_prmNumSecndry > 0)
            {
                double mcMainE   = 0.;
                double mcSecondE = 0.;

                // Identify the "primary" and "secondary" tracks in the conversion
                // Primary = highest energy particle
                // Secondary = lower energy of pair
                for(int idx = 0; idx < m_prmNumSecndry; idx++)
                {
                    const Event::McParticle* mcPart1 = m_mcEvent->getSecondary(idx);

                    // We are only interested in the tracks resulting from gamma conversion
                    // Note, however, that if process was purely compton, etc., then we will pass this
                    if (classifyBits & Event::McEventStructure::TRKCONVERT && mcPart1->getProcess() != "conv") 
                        continue;

                    // If we have already found a candidate primary, check against secondary
                    if (mcMain)
                    {
                        mcSecond  = mcPart1;
                        mcSecondE = mcPart1->initialFourMomentum().e();

                        // Primary defined to be highest energy of pair
                        if (mcMainE < mcSecondE)
                        {
                            mcSecond  = mcMain;
                            mcSecondE = mcMainE;
                            mcMain    = mcPart1;
                            mcMainE   = mcPart1->initialFourMomentum().e();
                        }
                    }
                    else
                    {
                        mcMain  = mcPart1;
                        mcMainE = mcPart1->initialFourMomentum().e();
                    }
                }
            }
        }

        // Make sure conversion is in a tracker tower
        if (curVolId[0] == 0 && curVolId[3] == 1 && curVolId.size() > 3)
        {
            m_prmDecLayer = 19;

            if (curVolId.size() > 6)
            {
                int trayNum = curVolId[4];
                int botTop  = curVolId[6];

                m_prmDecInCnv = botTop;
                m_prmDecLayer = trayNum;
            }
        }

        // Make sure we have found a primary particle
        if (mcMain)
        {
            // To get the sum of the angles
            Hep3Vector       mainVec    = m_mcEvent->getTrackDirection(mcMain);
            HepLorentzVector main4mom   = mcMain->initialFourMomentum();
            double           mainMom    = main4mom.vect().mag();
            HepLorentzVector main4cls(mainMom*mainVec,main4mom.e());

            // Fill the angle related info
            calcMcAngleInfo(mcMain, mcPart, m_scd1Ang2Gam, m_scd1Cls2Gam, m_scd1McTrkRms, m_scd1Type,
                            m_scd1Energy, m_scd1CosDirX, m_scd1CosDirY, m_scd1CosDirZ, m_scd1NHits, m_scd1NClstrs,
                            m_scd1NGaps, m_scd11stGapSz, m_scd1NHits2Gp);

            // Look here at track energy loss (if possible)
            calcMcEnergyInfo(mcMain, m_scd1ELastHit, m_scd1RadELoss, m_scd1DeltaRay,
                             m_scd1NBrems, m_scd1NDeltas, m_scd1NDeltaHt, m_scd1AveRange, m_scd1MaxRange);

            // Want to compress into the double scdType
            //short int* typeWords  = reinterpret_cast<short int*>(&m_scd1Type);

            //typeWords[0] = partType;
            //typeWords[1] = 0; //unused
            //typeWords[2] = m_mcEvent->getTrackHitLayer(mcMain,   0);
            //typeWords[3] = m_mcEvent->getTrackHitLayer(mcMain, 100);
            m_scd1FirstLyr = m_mcEvent->getTrackHitLayer(mcMain,   0);
            m_scd1LastLyr  = m_mcEvent->getTrackHitLayer(mcMain, 100);

            if (m_scd1RadELoss > m_scd1Energy - m_scd1ELastHit)
            {
                //int jj=0;
            }

            // Testing at this point
            const Event::TkrPatCand* patCand = m_mcTracks->getBestTkrPatCand(mcMain);
            if (patCand)
            {
#ifdef DEBUG
                int numHits = m_mcTracks->getNumMcHits(patCand,mcMain);
                int jjj=0;
#endif
            }

            // If there are two gamma conversion tracks then do this part
            if (mcSecond)
            {
                Hep3Vector       secondVec = m_mcEvent->getTrackDirection(mcSecond);
                HepLorentzVector scnd4mom  = mcSecond->initialFourMomentum();
                double           scndMom   = scnd4mom.vect().mag();
                HepLorentzVector scnd4cls(scndMom*secondVec,scnd4mom.e());
                
                calcMcAngleInfo(mcSecond, mcPart, m_scd2Ang2Gam, m_scd2Cls2Gam, m_scd2McTrkRms, m_scd2Type,
                                m_scd2Energy, m_scd2CosDirX, m_scd2CosDirY, m_scd2CosDirZ, m_scd2NHits, m_scd2NClstrs,
                                m_scd2NGaps, m_scd21stGapSz, m_scd2NHits2Gp);
                
                calcMcEnergyInfo(mcSecond, m_scd2ELastHit, m_scd2RadELoss, m_scd2DeltaRay,
                                 m_scd2NBrems, m_scd2NDeltas, m_scd2NDeltaHt, m_scd2AveRange, m_scd2MaxRange);
                
                m_scd2FirstLyr = m_mcEvent->getTrackHitLayer(mcSecond,   0);
                m_scd2LastLyr  = m_mcEvent->getTrackHitLayer(mcSecond, 100);

                // Check information on shared hits with other secondary tracks
                unsigned int sharedInfo = m_mcEvent->getSharedHitInfo(mcMain,mcSecond);
                //int          nShared    = sharedInfo & 0x7;

                m_prmTrkPattern = sharedInfo >> 4;

                // Update the MC direction from the four vector
                main4mom += scnd4mom;
                main4cls += scnd4cls;
            }

            // Look at reproducing the original gamma direction
            double cosMcGamAng  = main4mom.vect().unit().dot(primary4mom.vect().unit());
            double cosClsGamAng = main4cls.vect().unit().dot(primary4mom.vect().unit());

            if (cosMcGamAng > 1.0 ) cosMcGamAng = 1.0;
            if (cosClsGamAng > 1.0) cosClsGamAng = 1.0;

            m_prmMcAngle = acos(cosMcGamAng);
            m_prmClsAngle = acos(cosClsGamAng);

            if (m_prmMcAngle > 0.0001) 
            {
                //int checkithere = 0;
            }
        }
        else
        {
            //int ijk = 0;
        }
    }
    
    log << MSG::DEBUG << " returning. " << endreq;

    return sc;
}


void McAnalValsTool::calcMcAngleInfo(const Event::McParticle* mcPart, const Event::McParticle* mcGamma,
                                     double& ang2Gam, double& cls2Gam, double& mcTrkRms, double& prtType,
                                     double& prtEnergy, double& prtCosDirX, double& prtCosDirY, double& prtCosDirZ,
                                     double& prtNHits, double& prtNClstrs, double& prtNGaps, double& prt1stGapSz,
                                     double& prtNHits2Gp)
{
    Hep3Vector       mainVec    = m_mcEvent->getTrackDirection(mcPart);
    HepLorentzVector main4mom   = mcPart->initialFourMomentum();
    double           mainMom    = main4mom.vect().mag();
    HepLorentzVector main4cls(mainMom*mainVec,main4mom.e());
    double           cosAng2Gam = main4mom.vect().unit().dot(mcGamma->initialFourMomentum().vect().unit());

    if (ang2Gam < -1.) 
    {
        ang2Gam = -1.;
    }
    if (ang2Gam >  1.)
    {
        ang2Gam = 1.;
    }
    ang2Gam     = acos(cosAng2Gam);

    cosAng2Gam  = main4cls.vect().unit().dot(mcPart->initialFourMomentum().vect().unit());
    cls2Gam     = acos(cosAng2Gam);

    mcTrkRms    = m_mcEvent->getTrackStraightness(mcPart);

    prtType     = mcPart->particleProperty();
    prtEnergy   = main4mom.e();
    prtCosDirX  = main4mom.vect().unit().x();
    prtCosDirY  = main4mom.vect().unit().y();
    prtCosDirZ  = main4mom.vect().unit().z();

    prtNHits    = (m_mcEvent->getMcPartTrack(mcPart)).size();
    prtNClstrs  = m_mcEvent->getNumClusterHits(mcPart);
    prtNGaps    = m_mcEvent->getNumGaps(mcPart);

    //If track gaps, find size and check to see if first one shortens the track
    if (prtNGaps > 0)
    {
        prt1stGapSz = m_mcEvent->getGapSize(mcPart, 0);
        prtNHits2Gp = m_mcEvent->getGapStartHitNo(mcPart, 0) - 1;
    }

    return;
}

void McAnalValsTool::calcMcEnergyInfo(const Event::McParticle* mcPart,
                                      double& lastHitE, double& radLossE, double& dltaRayE, double& nBrems,
                                      double& nDeltas, double& nDeltaHt, double& aveRange, double& maxRange)
{
    int nBremsInt,nDeltasInt,nDeltaHtInt;

    lastHitE = m_mcEvent->getTrackELastHit(mcPart);
    radLossE = m_mcEvent->getTrackBremELoss(mcPart, nBremsInt);
    dltaRayE = m_mcEvent->getTrackDeltaELoss(mcPart, nDeltasInt, nDeltaHtInt);
    nDeltas  = m_mcEvent->getTrackDeltaRange(mcPart, aveRange, maxRange);

    nBrems   = nBremsInt;
    nDeltas  = nDeltasInt;
    nDeltaHt = nDeltaHtInt;

    return;
}

