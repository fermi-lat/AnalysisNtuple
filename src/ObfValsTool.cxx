/** @file ObfValsTool.cxx
@brief Calculates the Onboard Filter variables
@authors

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ObfValsTool.cxx,v 1.12.60.1 2008/04/28 21:52:58 usher Exp $
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

    // "Old" FilterStatus related variables
    double m_statusHi;
    double m_statusLo;
    double m_separation;
    double m_filtxdir,m_filtydir ,m_filtzdir;
    float  m_energy;
    double m_slopeYZ;
    double m_slopeXZ;
    int    m_xHits;
    int    m_yHits;

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
    int    m_hfcStatus;
    int    m_mipStatus;
    int    m_dfcStatus;

    int    m_gamEnergy;

    int    m_statusBytes;

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

Here is a list of the bits in the Gamma Filter word, as of March 2007:

@verbatim
    Name                          Bit  Explanation

    DFC_V_STATUS_ACD                   0  ACD was analyzed                
    DFC_V_STATUS_DIR                   1     DIR was decoded                 
    DFC_V_STATUS_ATF                   2     ACD/TKR veto was analyzed       
    DFC_V_STATUS_CAL1                  3  CAL was analyzed, phase 1       
    DFC_V_STATUS_TKR                   4  TKR finding was done            
    DFC_V_STATUS_ACD_TOP               5  ACD top  tile struck            
    DFC_V_STATUS_ACD_SIDE              6  ACD side tile struck            
    DFC_V_STATUS_ACD_SIDE_FILTER       7     ACD      filter tile struck     
    DFC_V_STATUS_TKR_POSSIBLE          8     Possible track                  
    DFC_V_STATUS_TKR_TRIGGER           9     Have a 3-in-a-row trigger       
    DFC_V_STATUS_CAL_LO               10  Cal Lo Trigger                  
    DFC_V_STATUS_CAL_HI               11  Cal Hi Trigger                  
    DFC_V_STATUS_TKR_EQ_1             12  Exactly 1 track                 
    DFC_V_STATUS_TKR_GE_2             13     Greater or equal 2 tracks       
    DFC_V_STATUS_TKR_THROTTLE         14     Throttle bit set                

    DFC_V_STATUS_TKR_LT_2_ELO         15  Low energy, no 2 track evidence   
    DFC_V_STATUS_TKR_SKIRT            16  Event into the skirt region     
    DFC_V_STATUS_TKR_EQ_0             17  No tracks                       
    DFC_V_STATUS_TKR_ROW2             18  Track Row 2 match               
    DFC_V_STATUS_TKR_ROW01            19  Track Row 0 or 1 match          
    DFC_V_STATUS_TKR_TOP              20  Track Top match                 
    DFC_V_STATUS_ZBOTTOM              21  No tracks into CAL with energy  
    DFC_V_STATUS_EL0_ETOT_90          22  E layer 0/ETOT > .90            
    DFC_V_STATUS_EL0_ETOT_01          23  E layer 0/ETOT < .01            
    DFC_V_STATUS_SIDE                 24     Event has a side face veto      
    DFC_V_STATUS_TOP                  25     Event has a top  face veto      
    DFC_V_STATUS_SPLASH_1             26  Event has a splash veto         
    DFC_V_STATUS_E350_FILTER_TILE     27  Event <350Mev  + filter tiles   
    DFC_V_STATUS_E0_TILE              28  Event 0 energy + tile hit       
    DFC_V_STATUS_SPLASH_0             29  Event has a splash veto         
    DFC_V_STATUS_NOCALLO_FILTER_TILE  30  No CAL LO trigger + filter tile 
    DFC_V_STATUS_VETOED               31  Any veto                        
@endverbatim

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td>  ObfGamStatus
<td>D<td>    The 32 bit status word output from the Gamma Filter 
<tr><td>  ObfHfcStatus
<td>D<td>    The 32 bit status word output from the Heavy Ion Filter 
<tr><td>  ObfMipStatus
<td>D<td>    The 32 bit status word output from the Minimum Ionizing Filter 
<tr><td>  ObfDfcStatus
<td>D<td>    The 32 bit status word output from the Diagnostic Filter 
<tr><td>  ObfStatusBytes
<td>D<td>    A 32 bit word which packs in the "status byte" output from each filter
<td>D<td>    Each filter contributes a 4 bit status word which describes "why" a given
<td>D<td>    filter sent an event into the data stream. Each 4 bit field is shifted into
<td>D<td>    a location by a factor given in the TDS Status class for that filter (see
<td>D<td>    the code below). 
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
<tr><td>  FilterStatus_HI
<td>D<td>    (deprecated) bits 15-31 of the filter status word (17 bits) 
<tr><td>  FilterStatus_LO
<td>D<td>    (deprecated) bits  0-14 of the filter status word (15 bits)
<tr><td>  FilterAngSep   
<td>D<td>    (deprecated) Filter status separation
<tr><td>  FilterEnergy
<td>F<td>    (deprecated) Energy as determined by onboard alg
<tr><td>  Filter
<td>I<td>    (deprecated) number of hits on best track XZ slope
<tr><td>  FilterYhits
<td>I<td>    (deprecated) number of hits on best track YZ slope
<tr><td>  FilterXZslope
<td>D<td>    (deprecated) XZ slope of the track
<tr><td>  FilterYZslope
<td>D<td>    (deprecated) YZ slope of the best track
<tr><td>  Filter[X/Y/Z]Dir
<td>D<td>    (deprecated) Direction cosines of the best track
</table> 

*/    


StatusCode ObfValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // Deprecated "old" FilterStatus variables
    addItem("FilterStatus_HI", &m_statusHi);
    addItem("FilterStatus_LO", &m_statusLo );
    addItem("FilterAngSep",    &m_separation );
    addItem("FilterEnergy",    &m_energy );
    addItem("FilterXhits",     &m_xHits );
    addItem("FilterYhits",     &m_yHits );
    addItem("FilterSlopeYZ",   &m_slopeYZ );
    addItem("FilterSlopeXZ",   &m_slopeXZ );
    addItem("FilterXDir",      &m_filtxdir );
    addItem("FilterYDir",      &m_filtydir );
    addItem("FilterZDir",      &m_filtzdir );

    // Status information from filters
    addItem("ObfGamStatus",    &m_gamStatus);
    addItem("ObfHfcStatus",    &m_hfcStatus);
    addItem("ObfMipStatus",    &m_mipStatus);
    addItem("ObfDfcStatus",    &m_dfcStatus);
    addItem("ObfStatusBytes",  &m_statusBytes);
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
            m_gamEnergy    = dynamic_cast<const OnboardFilterTds::ObfGammaStatus*>(obfResult)->getEnergy();
            m_statusBytes |= (obfResult->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::GammaFilter);

            filterStatusHi = m_gamStatus >> GFC_STATUS_V_TKR_LT_2_ELO;
            filterStatusLo = m_gamStatus && 0x7FFF;
        }
        else m_gamStatus = -1;

        // HFC Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::HFCFilter))
        {
            m_hfcStatus    = obfResult->getStatusWord();
            m_statusBytes |= (obfResult->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::HFCFilter);
        }
        else m_hfcStatus = -1;

        // MIP Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::MipFilter))
        {
            m_mipStatus    = obfResult->getStatusWord();
            m_statusBytes |= (obfResult->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::MipFilter);
        }
        else m_mipStatus = -1;

        // DFC Filter results
        if (obfResult   = 
                obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::DFCFilter))
        {
            m_dfcStatus    = obfResult ? obfResult->getStatusWord() : -1;
            m_statusBytes |= (obfResult->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::DFCFilter);
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

    // Copy deprecated variables from new stuff
    m_statusHi   = filterStatusHi;
    m_statusLo   = filterStatusLo;
    m_energy     = m_gamEnergy;
    m_xHits      = m_nGrbHitsX;
    m_yHits      = m_nGrbHitsY;
    m_filtxdir   = m_grbXdir;
    m_filtydir   = m_grbYdir;
    m_filtzdir   = m_grbZdir;
    m_slopeYZ    = m_grbSlpX;
    m_slopeXZ    = m_grbSlpY;
    m_separation = m_grbAngSep;

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

