/** @file ObfValsTool.cxx
@brief Calculates the Onboard Filter variables
@authors

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ObfValsTool.cxx,v 1.10 2007/12/29 18:48:29 usher Exp $
*/

// Include files

#include "ValBase.h"

//#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
//#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/MonteCarlo/McParticle.h"

#include "CLHEP/Vector/Rotation.h"
#include "geometry/Vector.h"

#include "OnboardFilterTds/FilterStatus.h"
#include "OnboardFilterTds/ObfFilterStatus.h"
#include "OnboardFilterTds/FilterAlgTds.h"

// to get current position
//flux
//#include "FluxSvc/IFluxSvc.h"
//#include "astro/GPS.h"

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

    double m_statusHi;
    double m_statusLo;
    double m_separation;
    double m_filterAlgStatus;
    double m_filtxdir,m_filtydir ,m_filtzdir;
    double m_grbXdir;
    double m_grbYdir;
    double m_grbZdir;
    double m_grbAngSep;

    //float OBF_ra, OBF_dec; // set by astro::GPS 
    //float OBF_glon, OBF_glat;

    float  m_energy;
    double m_slopeYZ;
    double m_slopeXZ;
    int    m_xHits;
    int    m_yHits;

    double m_grbSlpX;
    double m_grbSlpY;
    int    m_nGrbHitsX;
    int    m_nGrbHitsY;

    int    m_gamStatus;
    int    m_hfcStatus;
    int    m_mipStatus;
    int    m_dfcStatus;

    int    m_warnNoFilterStatus;

    //IFluxSvc*   m_fluxSvc;
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

Here is a list of the bits in the filter word, as of March 2007:

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
<tr><td>  FilterAlgStatus
<td>D<td>    Status generated by FilterAlg 
<tr><td>  FilterStatus_HI
<td>D<td>    bits 15-31 of the filter status word (17 bits) 
<tr><td>  FilterStatus_LO
<td>D<td>    bits  0-14 of the filter status word (15 bits)
<tr><td>  FilterAngSep   
<td>D<td>    Filter status separation
<tr><td>  FilterEnergy
<td>F<td>    Energy as determined by onboard alg
<tr><td>  Filter
<td>I<td>    number of hits on best track XZ slope
<tr><td>  FilterYhits
<td>I<td>    number of hits on best track YZ slope
<tr><td>  FilterXZslope
<td>D<td>    XZ slope of the track
<tr><td>  FilterYZslope
<td>D<td>    YZ slope of the best track
<tr><td>  Filter[X/Y/Z]Dir
<td>D<td>    Direction cosines of the best track
<tr><td>  Obf[Gam/Hfc/Mip]Status
<td>I<td>    Status bits for the Gamma, HFC and MIP filters
</table> 

*/    


StatusCode ObfValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the services    
    //if ( service("FluxSvc", m_fluxSvc).isFailure() ){
    //    log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
    //    return StatusCode::FAILURE;
    //}

    addItem("FilterStatus_HI", &m_statusHi);
    addItem("FilterStatus_LO", &m_statusLo );
    addItem("FilterAlgStatus", &m_filterAlgStatus );
    addItem("FilterAngSep",    &m_separation );
    addItem("FilterEnergy",    &m_energy );
    addItem("FilterXhits",     &m_xHits );
    addItem("FilterYhits",     &m_yHits );
    addItem("FilterSlopeYZ",   &m_slopeYZ );
    addItem("FilterSlopeXZ",   &m_slopeXZ );
    addItem("FilterXDir",      &m_filtxdir );
    addItem("FilterYDir",      &m_filtydir );
    addItem("FilterZDir",      &m_filtzdir );

    addItem("GrbXHits",        &m_nGrbHitsX);
    addItem("GrbYHits",        &m_nGrbHitsY);
    addItem("GrbSlpX",         &m_grbSlpX);
    addItem("GrbSlpY",         &m_grbSlpY);
    addItem("GrbXDir",         &m_grbXdir);
    addItem("GrbYDir",         &m_grbYdir);
    addItem("GrbZDir",         &m_grbZdir);
    addItem("GrbMcAngSep",     &m_grbAngSep);

    //addItem("FilterxRa",        &OBF_ra);
    //addItem("FilterxDec",       &OBF_dec); 
    //addItem("FilterxL",         &OBF_glon);
    //addItem("FilterxB",         &OBF_glat);

    addItem("ObfGamStatus",    &m_gamStatus);
    addItem("ObfHfcStatus",    &m_hfcStatus);
    addItem("ObfMipStatus",    &m_mipStatus);
    addItem("ObfDfcStatus",    &m_dfcStatus);

    zeroVals();

    m_warnNoFilterStatus = 0;

    return sc;
}

StatusCode ObfValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Old school output
    SmartDataPtr<OnboardFilterTds::FilterStatus> 
        filterStatus(m_pEventSvc, "/Event/Filter/FilterStatus");
    if( filterStatus )
    {
        m_statusHi   = filterStatus->getHigh();
        m_statusLo   = filterStatus->getLow();
        m_separation = filterStatus->getSeparation();
        m_energy     = filterStatus->getCalEnergy();
        m_xHits      = 0;
        m_yHits      = 0;
        m_filtxdir   = m_filtydir = m_filtzdir = 0;

        double slopeXZ  = 0.0;
        double slopeYZ  = 0.0;
        double intXZ    = 0.0;
        double intYZ    = 0.0;

        filterStatus->getBestTrack(m_xHits,m_yHits,slopeXZ,slopeYZ,intXZ,intYZ, m_nGrbHitsX, m_nGrbHitsY, m_grbSlpX, m_grbSlpY);

        // Look up the McParticle collection in case its there... 
        SmartDataPtr<Event::McParticleCol> mcParticleCol(m_pEventSvc, EventModel::MC::McParticleCol);

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
        
        if(m_xHits > 0 && m_yHits > 0)
        {
            float alpha = atan2(slopeYZ,slopeXZ);
            if(alpha < 0) 
            {
                alpha = alpha+2.0*3.1415;
            }
            float m_slope = sqrt(pow(slopeXZ,2) + pow(slopeYZ,2));
            float beta    = atan(m_slope);

            m_filtxdir = cos(alpha)*sin(beta);
            m_filtydir = sin(alpha)*sin(beta);
            m_filtzdir = cos(beta);

            // code below moved to ObfCoordsAlg to accomodate Interleave
            /*
            Vector filtDir(m_filtxdir, m_filtydir, m_filtzdir);

            // convert to (ra, dec)

            // The GPS singleton has current time and orientation
            static astro::GPS* gps = m_fluxSvc->GPSinstance();
            double time = gps->time();

            CLHEP::HepRotation R ( gps->transformToGlast(time, astro::GPS::CELESTIAL) );

            // no minus sign... filter tracks already point up!

            astro::SkyDir skydir( R.inverse() * filtDir );
            OBF_ra   = skydir.ra();
            OBF_dec  = skydir.dec();
            OBF_glon = skydir.l();
            OBF_glat = skydir.b();
            */

            // If running Monte Carlo, determine FilterAngSep here
            if (mcParticleCol)
            {
                // We only care about incident particle direction here, so can use the first particle
                // in the McParticleCol
                Event::McParticle* mcPart = *(mcParticleCol->begin());

                CLHEP::HepLorentzVector Mc_p0 = mcPart->initialFourMomentum();
                Vector                  Mc_t0 = Vector(Mc_p0.x(),Mc_p0.y(), Mc_p0.z()).unit();
                Vector                  filtDir(-m_filtxdir, -m_filtydir, -m_filtzdir);

                double cosTheta = filtDir.dot(Mc_t0);

                m_separation = acos(cosTheta);
            }
        }

        m_slopeYZ = slopeYZ;
        m_slopeXZ = slopeXZ;
    }else {
        m_statusHi = m_statusLo = 0;

        m_warnNoFilterStatus++;
        if (   m_warnNoFilterStatus <= 10 ) {
            log << MSG::WARNING << "FilterStatus not found" ;
            if ( m_warnNoFilterStatus == 10 ) {
                log << " -- Further WARNINGs on missing FilterStatus are suppressed"; }
            log  << endreq;
        }
    }

    // Beyond old school
    SmartDataPtr<FilterAlgTds::FilterAlgData> 
        filterAlgStatus(m_pEventSvc,"/Event/Filter/FilterAlgData");
    if(filterAlgStatus){
        m_filterAlgStatus=(double)filterAlgStatus->getVetoWord();
    }

    // ultra modern method
    SmartDataPtr<OnboardFilterTds::ObfFilterStatus> 
        obfStatus(m_pEventSvc, "/Event/Filter/ObfFilterStatus");

    if (obfStatus)
    {
        // Pointer to our retrieved objects
        const OnboardFilterTds::IObfStatus* obfResult = 0;

        // Start with Gamma Filter
        obfResult   = 
            obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::GammaFilter);
        m_gamStatus = obfResult ? obfResult->getStatus32() : -1;

        // HFC Filter results
        obfResult   = 
            obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::HFCFilter);
        m_hfcStatus = obfResult ? obfResult->getStatus32() : -1;

        // MIP Filter results
        obfResult   = 
            obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::MipFilter);
        m_mipStatus = obfResult ? obfResult->getStatus32() : -1;

        // DFC Filter results
        obfResult   = 
            obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::DFCFilter);
        m_dfcStatus = obfResult ? obfResult->getStatus32() : -1;
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

