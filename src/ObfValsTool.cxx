/** @file ObfValsTool.cxx
@brief Calculates the Onboard Filter variables
@authors

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ObfValsTool.cxx,v 1.13 2008/05/21 16:58:52 usher Exp $
*/

// Include files

#include "ValBase.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/MonteCarlo/McParticle.h"

#include "CLHEP/Vector/Rotation.h"
#include "geometry/Vector.h"

#include "OnboardFilterTds/ObfFilterStatus.h"
#include "OnboardFilterTds/ObfFilterTrack.h"

/*! @class ObfValsTool
@brief calculates Obf values

@authors 
*/

class ObfValsTool : public ValBase
{
public:

    ObfValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~ObfValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

    StatusCode finalize();

private:

    // "Best" track variables
    int    m_nGrbHitsX;
    int    m_nGrbHitsY;
    float  m_grbSlpX;
    float  m_grbSlpY;
    float  m_grbIntX;
    float  m_grbIntY;
    float  m_grbZ;
    float  m_grbXdir;
    float  m_grbYdir;
    float  m_grbZdir;
    float  m_grbAngSep;

    // Filter status information
    int    m_gamStatus;
    int    m_gamState;
    int    m_hfcStatus;
    int    m_hfcState;
    int    m_mipStatus;
    int    m_mipState;
    int    m_dfcStatus;
    int    m_dfcState;

    int    m_gamStage;
    float  m_gamEnergy;

    int    m_warnNoFilterStatus;
};

// Static factory for instantiation of algtool objects
static ToolFactory<ObfValsTool> s_factory;
const IToolFactory& ObfValsToolFactory = s_factory;

// Standard Constructor
ObfValsTool::ObfValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

/** @page anatup_vars 
@section obfvalstool ObfValsTool Variables

The bit definitions for the output status words are contained in the following files 
located in the Glast obf external package (starting at obf/B1-0-8/PHY):
1) Gamma Filter - EFC/GFC_status.h
2) HIP Filter - XFC/HFC_status.h
3) MIP Filter - XFC/MFC_status.h
4) DGN Filter - XFC/DFC_status.h

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td>  ObfGamStatus
<td>D<td>    The 32 bit status word output from the Gamma Filter 
<tr><td>  ObfGamState
<td>D<td>    The 4 bit state word summarizes veto and prescaler states 
<tr><td>  ObfHfcStatus
<td>D<td>    The 32 bit status word output from the Heavy Ion Filter 
<tr><td>  ObfHfcState
<td>D<td>    The 4 bit state word summarizes veto and prescaler states 
<tr><td>  ObfMipStatus
<td>D<td>    The 32 bit status word output from the Minimum Ionizing Filter 
<tr><td>  ObfMipState
<td>D<td>    The 4 bit state word summarizes veto and prescaler states 
<tr><td>  ObfDfcStatus
<td>D<td>    The 32 bit status word output from the Diagnostic Filter 
<tr><td>  ObfDfcState
<td>D<td>    The 4 bit state word summarizes veto and prescaler states 
<tr><td>  ObfGamStage
<tr>D<td>    The "stage" reached by the gamma filter (see the GFC_STAGE definition
             in EFC/GFC_status.h)
<tr><td>  ObfGamEnergy
<td>D<td>    The energy in the calorimeter seen by the Gamma Filter
<tr><td>  GrbXHits,GrbYHits
<td>D<td>    Number of hits in X (or Y) projection of "best" track
<tr><td>  GrbSlpX,GrbSlpY
<td>D<td>    Slope of "best" track in X or Y projection
<tr><td>  GrbIntX,GrbIntY
<td>D<td>    Intercept in X-Z or Y-Z projection of "best" track
<tr><td>  GrbZ
<td>D<td>    Z coordinate of "best" track
<tr><td>  GrbXDir,GrbYDir,GrbZDir
<td>D<td>    Direction cosines of "best" track
<tr><td>  GrbMcAngSep
<td>D<td>    If MC run, angle to "true" incident particle of "best" track
</table> 

*/    


StatusCode ObfValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // Status information from filters
    addItem("ObfGamStatus",    &m_gamStatus);
    addItem("ObfGamState",     &m_gamState);
    addItem("ObfHipStatus",    &m_hfcStatus);
    addItem("ObfHipState",     &m_hfcState);
    addItem("ObfMipStatus",    &m_mipStatus);
    addItem("ObfMipState",     &m_mipState);
    addItem("ObfDgnStatus",    &m_dfcStatus);
    addItem("ObfDgnState",     &m_dfcState);

    addItem("ObfGamStage",     &m_gamStage);
    addItem("ObfGamEnergy",    &m_gamEnergy);

    // "Best" track information
    addItem("GrbXHits",        &m_nGrbHitsX);
    addItem("GrbYHits",        &m_nGrbHitsY);
    addItem("GrbSlpX",         &m_grbSlpX);
    addItem("GrbSlpY",         &m_grbSlpY);
    addItem("GrbIntX",         &m_grbIntX);
    addItem("GrbIntY",         &m_grbIntY);
    addItem("GrbZ",            &m_grbZ);
    addItem("GrbXDir",         &m_grbXdir);
    addItem("GrbYDir",         &m_grbYdir);
    addItem("GrbZDir",         &m_grbZdir);
    addItem("GrbMcAngSep",     &m_grbAngSep);

    zeroVals();

    m_warnNoFilterStatus = 0;

    return sc;
}

StatusCode ObfValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // We'll use these when copying to the deprecated variables
    unsigned int filterStatusHi = 0;
    unsigned int filterStatusLo = 0;

    // Start with the status information returned from the filters
    // Look up the TDS class containing this information
    SmartDataPtr<OnboardFilterTds::ObfFilterStatus> 
        obfStatus(m_pEventSvc, "/Event/Filter/ObfFilterStatus");

    // If it exists then fill filter status info
    if (obfStatus)
    {
        // Pointer to our retrieved objects
        const OnboardFilterTds::IObfStatus* obfResult = 0;

        // Start with Gamma Filter
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::GammaFilter))
        {
            m_gamStatus    = obfResult->getStatusWord();
            m_gamState     = obfResult->getFiltersb() >> 4;
            m_gamStage     = dynamic_cast<const OnboardFilterTds::ObfGammaStatus*>(obfResult)->getStage();
            m_gamEnergy    = dynamic_cast<const OnboardFilterTds::ObfGammaStatus*>(obfResult)->getEnergy() / 4.;
        }
        else m_gamStatus = -1;

        // HFC Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::HIPFilter))
        {
            m_hfcStatus    = obfResult->getStatusWord();
            m_hfcState     = obfResult->getFiltersb() >> 4;
        }
        else m_hfcStatus = -1;

        // MIP Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::MIPFilter))
        {
            m_mipStatus    = obfResult->getStatusWord();
            m_mipState     = obfResult->getFiltersb() >> 4;
        }
        else m_mipStatus = -1;

        // DFC Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::DGNFilter))
        {
            m_dfcStatus    = obfResult ? obfResult->getStatusWord() : -1;
            m_dfcState     = obfResult->getFiltersb() >> 4;
        }
        else m_dfcStatus = -1;
    }
    else 
    {
        m_warnNoFilterStatus++;
        if (   m_warnNoFilterStatus <= 10 ) 
        {
            log << MSG::WARNING << "ObfFilterStatus not found" ;
            if ( m_warnNoFilterStatus == 10 ) {
                log << " -- Further WARNINGs on missing FilterStatus are suppressed"; }
            log  << endreq;
        }
    }

    // Get the summary tracking information on the "best" track
    SmartDataPtr<OnboardFilterTds::ObfFilterTrack> 
        filterTrack(m_pEventSvc, "/Event/Filter/ObfFilterTrack");

    if( filterTrack )
    {
        m_nGrbHitsX = filterTrack->get_nXhits();
        m_nGrbHitsY = filterTrack->get_nYhits();
        m_grbSlpX   = filterTrack->get_slpXZ();
        m_grbSlpY   = filterTrack->get_slpYZ();
        m_grbIntX   = filterTrack->get_xInt();
        m_grbIntY   = filterTrack->get_yInt();
        m_grbZ      = filterTrack->get_z();

        // Calculate direction cosines for GRB track first
        if (m_nGrbHitsX > 0 && m_nGrbHitsY > 0)
        {
            double alpha = atan2(m_grbSlpY, m_grbSlpX);

            if (alpha < 0.) alpha += 2. * M_PI;

            double slope = sqrt(pow(m_grbSlpX,2) + pow(m_grbSlpY,2));
            double beta  = atan(slope);

            m_grbXdir = cos(alpha)*sin(beta);
            m_grbYdir = sin(alpha)*sin(beta);
            m_grbZdir = cos(beta);

            // Look up the McParticle collection in case its there... 
            SmartDataPtr<Event::McParticleCol> mcParticleCol(m_pEventSvc, EventModel::MC::McParticleCol);

            // If running Monte Carlo, determine FilterAngSep here
            if (mcParticleCol)
            {
                // We only care about incident particle direction here, so can use the first particle
                // in the McParticleCol
                Event::McParticle* mcPart = *(mcParticleCol->begin());

                CLHEP::HepLorentzVector Mc_p0 = mcPart->initialFourMomentum();
                Vector                  Mc_t0 = Vector(Mc_p0.x(),Mc_p0.y(), Mc_p0.z()).unit();
                Vector                  filtDir(-m_grbXdir, -m_grbYdir, -m_grbZdir);

                double cosTheta = filtDir.dot(Mc_t0);

                m_grbAngSep = acos(cosTheta);
            }
        }
    }

    return sc;
}

StatusCode ObfValsTool::finalize() {

    MsgStream log(msgSvc(), name());
    log << MSG::INFO ;
    log << endreq;
    if(m_warnNoFilterStatus>0)
        log << MSG::INFO 
        << "Number of warnings (FilterStatus not found): "
        << m_warnNoFilterStatus << endreq;

    //setFinalized(); //  prevent being called again

    return StatusCode::SUCCESS;
}

