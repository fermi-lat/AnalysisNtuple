/** @file ObfValsTool.cxx
@brief Calculates the Onboard Filter variables
@authors

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ObfValsTool.cxx,v 1.23 2008/08/01 04:13:41 heather Exp $
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
#include "LdfEvent/LsfMetaEvent.h"

/*! @class ObfValsTool
@brief calculates Obf values

@authors 
*/

namespace {
    // save a bit of typing!
    int _invalid = enums::Lsf::INVALID;
}

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

    void zeroVals();

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
    int    m_hipStatus;
    int    m_hipState;
    int    m_mipStatus;
    int    m_mipState;
    int    m_dgnStatus;
    int    m_dgnState;

    int    m_gamStage;
    float  m_gamEnergy;

    int    m_warnNoFilterStatus;

// Fsw filter info
    int    m_fswGamStatus;
    int    m_fswGamState;
    int    m_fswGamPrescaleIndex;
    int    m_fswGamPrescaleFactor;
    float  m_fswGamEnergy;
    int    m_fswGamStage;
    int    m_fswGamVersion;

    int    m_fswHipStatus;
    int    m_fswHipState;
    int    m_fswMipStatus;
    int    m_fswMipState;

    int    m_fswDgnStatus;
    int    m_fswDgnState;
    int    m_fswDgnPrescaleIndex;
    int    m_fswDgnPrescaleFactor;

    void copyObfToFsw();
  
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
-# Gamma Filter - EFC/GFC_status.h
-# HIP Filter - XFC/HFC_status.h
-# MIP Filter - XFC/MFC_status.h
-# DGN Filter - XFC/DFC_status.h

The ObfXXX and GrbXXX variables are calculated in the ground-software
simulation of the onboard flight software.

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td>  ObfGamStatus, same for Hip, Mip, Dgn</td>
<td>I<td>    The 32-bit status word output from the Gamma, Heavy-Ion, 
             Minimum-Ionizing and Diagnostic Filters, resp, as calculated
             by the OnboardFilter simulator. </td>
<tr><td>  ObfGamState, same for Hip, Mip, Dgn</td>
<td>I<td>    The 4-bit state word summarizes veto and prescaler states, as
             calculated by the OnboardFilter simulator.
<tr><td>  ObfGamStage
<td>I<td>    The "stage" reached by the gamma filter (see the GFC_STAGE definition
             in EFC/GFC_status.h)
<tr><td>  ObfGamEnergy
<td>F<td>    The energy in the calorimeter seen by the Gamma Filter
<tr><td>  Grb[X/Y]Hits
<td>I<td>    Number of hits in X (or Y) projection of "best" track
<tr><td>  GrbSlp[X/Y]
<td>F<td>    Slope of "best" track in X or Y projection
<tr><td>  GrbInt[X/Y]
<td>F<td>    Intercept in X-Z or Y-Z projection of "best" track
<tr><td>  GrbZ
<td>F<td>    Z coordinate of "best" track
<tr><td>  Grb[X/Y/Z]Dir
<td>F<td>    Direction cosines of "best" track
<tr><td>  GrbMcAngSep
<td>F<td>    If MC run, angle to "true" incident particle of "best" track
</table>

The FswXXX variables are calculated onboard in flight software. For simulated data,
they will be copies of the corresponding ObfXXX variables.

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td>  FswGamStatus, same for Hip, Mip, Dgn
<td>I<td>    The 32-bit status word output from the Gamma, Heavy-Ion, 
             Minimum-Ionizing and Diagnostic Filters, resp.
<tr><td>  FswGamState, same for Hip, Mip, Dgn</td>
<td>I<td>    The 4-bit state word summarizes veto and prescaler states
<tr><td>  FswGamPrescaleIndex, same for Dgn</td>
<td>I<td>    The prescaler which expired for this event.   
             Remapped to the order they are applied, counting down from 33.<br>
                33 -> input prescaler.  Event is ignored by this filter.<br>
                32 -> output prescaler. Event is passed by this filter.<br>
                31-0 -> "line" prescalers.  Only checked when each bit is set.<br>
                -1 -> no prescaler expired.
<tr><td>  FswGamPrescaleFactor, same for Dgn </td>
<td>I<td>    The value to which the relevent prescaler counted 
             before it expired.  ie.  the number of similar events 
             that were tossed for each one that was kept, the "weight"
<tr><td>  FswGamEnergy</td>
<td>F<td>    See ObfGamEnergy
<tr><td>  FswGamStage</td>
<td>I<td>    See ObfGamStage
<tr><td>  FswGamVersion</td>
<td>I<td>   Version of the Gamma Filter 
</table> 

*/    

StatusCode ObfValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // Status information from filters
    addItem("ObfGamStatus",         &m_gamStatus);
    addItem("ObfGamState",          &m_gamState);
    addItem("ObfGamStage",          &m_gamStage);
    addItem("ObfGamEnergy",         &m_gamEnergy);

    addItem("ObfHipStatus",         &m_hipStatus);
    addItem("ObfHipState",          &m_hipState);

    addItem("ObfMipStatus",         &m_mipStatus);
    addItem("ObfMipState",          &m_mipState);

    addItem("ObfDgnStatus",         &m_dgnStatus);
    addItem("ObfDgnState",          &m_dgnState);

    // "Best" track information
    addItem("GrbXHits",             &m_nGrbHitsX);
    addItem("GrbYHits",             &m_nGrbHitsY);
    addItem("GrbSlpX",              &m_grbSlpX);
    addItem("GrbSlpY",              &m_grbSlpY);
    addItem("GrbIntX",              &m_grbIntX);
    addItem("GrbIntY",              &m_grbIntY);
    addItem("GrbZ",                 &m_grbZ);
    addItem("GrbXDir",              &m_grbXdir);
    addItem("GrbYDir",              &m_grbYdir);
    addItem("GrbZDir",              &m_grbZdir);
    addItem("GrbMcAngSep",          &m_grbAngSep);

    addItem("FswGamStatus",         &m_fswGamStatus);
    addItem("FswGamState",          &m_fswGamState);
    addItem("FswGamEnergy",         &m_fswGamEnergy);
    addItem("FswGamStage",          &m_fswGamStage);
    addItem("FswGamPrescaleIndex",  &m_fswGamPrescaleIndex);
    addItem("FswGamPrescaleFactor", &m_fswGamPrescaleFactor);
    addItem("FswGamVersion",        &m_fswGamVersion);

    addItem("FswHipStatus",         &m_fswHipStatus);
    addItem("FswHipState",          &m_fswHipState);
    addItem("FswMipStatus",         &m_fswMipStatus);
    addItem("FswMipState",          &m_fswMipState);

    addItem("FswDgnStatus",         &m_fswDgnStatus);
    addItem("FswDgnState",          &m_fswDgnState);
    addItem("FswDgnPrescaleIndex",  &m_fswDgnPrescaleIndex);
    addItem("FswDgnPrescaleFactor", &m_fswDgnPrescaleFactor);


    zeroVals();

    m_warnNoFilterStatus = 0;

    return sc;
}

void ObfValsTool::zeroVals() {
    
    // We don't want to zero these unless
    //  there's a digi TDS around to restore them from.
    // So we'll do the zeroing in calculate()

    return;
}

StatusCode ObfValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // We'll use these when copying to the deprecated variables
    unsigned int filterStatusHi = 0;
    unsigned int filterStatusLo = 0;


    // Obf variables
    SmartDataPtr<OnboardFilterTds::ObfFilterStatus> 
        obfStatus(m_pEventSvc, "/Event/Filter/ObfFilterStatus");
    // Fsw variables
    SmartDataPtr<LsfEvent::MetaEvent>  
        metaEventTds(m_pEventSvc, "/Event/MetaEvent");
    const lsfData::GammaHandler* gamma;
    const lsfData::DgnHandler*   dgn;
    const lsfData::HipHandler*   hip;       
    const lsfData::MipHandler*   mip;
    if(metaEventTds) {
        gamma = metaEventTds->gammaFilter();
        dgn   = metaEventTds->dgnFilter();
        hip   = metaEventTds->hipFilter();       
        mip   = metaEventTds->mipFilter();
    }

    // Grb variables
    SmartDataPtr<OnboardFilterTds::ObfFilterTrack> 
        filterTrack(m_pEventSvc, "/Event/Filter/ObfFilterTrack");

    // There is no digi Tds, 
    //  so copy Obf to Fsw if the latter is set to defaults
    bool isMetaEventData = gamma||dgn||mip||hip;
    if(!obfStatus && !isMetaEventData && !filterTrack) {
        copyObfToFsw();
        return sc;
    }
    
    // Start with the status information returned from the filters
    // Look up the TDS class containing this information

    if (obfStatus)
    {
        // If obfStatus exists then fill filter status info

        m_gamState = _invalid;
        m_dgnState = _invalid;
        m_hipState = _invalid;
        m_mipState = _invalid;
        m_gamStatus = _invalid;
        m_dgnStatus = _invalid;
        m_hipStatus = _invalid;
        m_mipStatus = _invalid;

        // Pointer to our retrieved objects
        const OnboardFilterTds::IObfStatus* obfResult = 0;

        // Start with Gamma Filter
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::GammaFilter))
        {
            m_gamStatus    = obfResult->getStatusWord();
            m_gamState     = obfResult->getState();
            m_gamStage     = 
                dynamic_cast<const OnboardFilterTds::ObfGammaStatus*>(obfResult)->getStage();
            m_gamEnergy    = 
                dynamic_cast<const OnboardFilterTds::ObfGammaStatus*>(obfResult)->getEnergy();
        }

        // HFC Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::HIPFilter))
        {
            m_hipStatus    = obfResult->getStatusWord();
            m_hipState     = obfResult->getState();
        }

        // MIP Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::MIPFilter))
        {
            m_mipStatus    = obfResult->getStatusWord();
            m_mipState     = obfResult->getState();
        }

        // DFC Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::DGNFilter))
        {
            m_dgnStatus    = obfResult ? obfResult->getStatusWord() : -1;
            m_dgnState     = obfResult->getState();
        }
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
    // Next, try to find the Fsw versions of the above

    if (metaEventTds) {
        if (gamma) {
            // rsd is Result Summary Data 
            if(gamma->rsd()) {
                m_fswGamStatus = gamma->rsd()->status();  
                m_fswGamVersion = gamma->version();
            } 
            m_fswGamState = gamma->state();  
            m_fswGamPrescaleFactor = gamma->prescaleFactor();    
            m_fswGamPrescaleIndex  = gamma->lpaHandler().prescaleIndex();
            m_fswGamEnergy         = static_cast<float>(gamma->rsd()->energyInLeus())/4.0;
            m_fswGamStage          = gamma->rsd()->stage();
        } else {
            m_fswGamPrescaleIndex = LSF_INVALID_UINT;
            m_fswGamPrescaleFactor = LSF_INVALID_UINT;
            m_fswGamState = _invalid;
            m_fswGamStatus = _invalid;
            m_fswGamEnergy = 0.0;
            m_fswGamStage = _invalid;
            m_fswGamVersion = _invalid;
        }

        if (dgn) {
            // rsd is Result Summary Data 
            if(dgn->rsd()) m_fswDgnStatus = dgn->rsd()->status();  
            m_fswDgnState = dgn->state();  
            m_fswDgnPrescaleFactor = dgn->prescaleFactor();  
            m_fswDgnPrescaleIndex  = dgn->lpaHandler().prescaleIndex(); 
        } else {
            m_fswDgnPrescaleIndex = LSF_INVALID_UINT;
            m_fswDgnPrescaleFactor = LSF_INVALID_UINT;
            m_fswDgnState = _invalid;
            m_fswDgnStatus = _invalid;
        }
 
        if (hip) {
            // rsd is Result Summary Data 
            if(hip->rsd()) m_fswHipStatus = hip->rsd()->status();  
            m_fswHipState = hip->state();  
            //m_fswHipPrescaleFactor = hip->prescaleFactor();  
            //m_fswHipPrescaleIndex  = hip->lpaHandler().prescaleIndex(); 
        } else {
            m_fswHipState = _invalid;
            m_fswHipStatus = _invalid;
        }

        if (mip) {
            // rsd is Result Summary Data 
            if(mip->rsd()) m_fswMipStatus = mip->rsd()->status();  
            m_fswMipState = mip->state();  
            //m_fswMipPrescaleFactor = mip->prescaleFactor();  
            //m_fswMipPrescaleIndex  = mip->lpaHandler().prescaleIndex(); 
        } else {
            m_fswMipState = _invalid;
            m_fswMipStatus = _invalid;
        }

        if(!gamma && !dgn && !mip && !hip) {
            log << MSG::DEBUG << "No MetaEvent" << endreq;
            // copy Obf to Fsw if appropriate
            copyObfToFsw();
        }
    }

    // Get the summary tracking information on the "best" track

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
            SmartDataPtr<Event::McParticleCol> 
                mcParticleCol(m_pEventSvc, EventModel::MC::McParticleCol);

            // If running Monte Carlo, determine FilterAngSep here
            m_grbAngSep = 0.0;
            if (mcParticleCol)
            {
                // We only care about incident particle direction here, 
                //   so can use the first particle in the McParticleCol
                Event::McParticle* mcPart = *(mcParticleCol->begin());

                CLHEP::HepLorentzVector Mc_p0 = mcPart->initialFourMomentum();

                Vector Mc_t0 = Vector(Mc_p0.x(),Mc_p0.y(), Mc_p0.z()).unit();
                Vector filtDir(-m_grbXdir, -m_grbYdir, -m_grbZdir);

                double cosTheta = filtDir.dot(Mc_t0);

                m_grbAngSep = acos(cosTheta);
            }
        } else {
            m_nGrbHitsX = 0;
            m_nGrbHitsY = 0;
            m_grbSlpX   = 0.0;
            m_grbSlpY   = 0.0;
            m_grbIntX   = 0.0;
            m_grbIntY   = 0.0;
            m_grbZ      = 0.0;

            m_grbXdir   = 0.0;
            m_grbYdir   = 0.0;
            m_grbZdir   = 0.0;
            m_grbAngSep = 0.0;
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

void ObfValsTool::copyObfToFsw()
{
    if(m_fswGamState == _invalid && m_fswDgnState == _invalid
        && m_fswHipState == _invalid && m_fswMipState == _invalid) 
    {

        m_fswGamStatus = m_gamStatus;
        m_fswGamState  = m_gamState;
        m_fswGamStage  = m_gamStage;
        m_fswGamEnergy = m_gamEnergy;
        m_fswMipStatus = m_mipStatus;
        m_fswMipState  = m_mipState;
        m_fswHipStatus = m_hipStatus;
        m_fswHipState  = m_hipState;
        m_fswDgnStatus = m_dgnStatus;
        m_fswDgnState  = m_dgnState;
    }
}

