/** @file OverlayValsTool.cxx
@brief declaration and definition of the class OverlayValsTool

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/OverlayValsTool.cxx,v 1.0 2008/10/01 07:02:16 lsrea Exp $

*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartDataLocator.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

// Event for access to time
#include "OverlayEvent/OverlayEventModel.h"
#include "OverlayEvent/EventOverlay.h"
#include "OverlayEvent/GemOverlay.h"
#include "OverlayEvent/PtOverlay.h"
#include "OverlayEvent/AcdOverlay.h"
#include "OverlayEvent/CalOverlay.h"
#include "OverlayEvent/TkrOverlay.h"

// stuff for exposure info (not really MC)
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/Exposure.h"

// to write a Tree with pointing info
#include "ntupleWriterSvc/INTupleWriterSvc.h"


#include "facilities/Util.h"
#include "facilities/Timestamp.h"

#include "astro/JulianDate.h"
#include "astro/GPS.h"

//
#include "AnalysisNtuple/PointingInfo.h"

namespace { // anonymous namespace for file-global
    astro::GPS* gps(0);  // pointer to relevant GPS entry
    const double R2D = 180./M_PI;
}

using namespace AnalysisNtuple;

/** 
* \class OverlayValsTool
*
* \brief This is an Algorithm designed to store pointing information in the tuple
* \author Toby Burnett
* 
*/


class OverlayValsTool  : public ValBase 
{
public:
    OverlayValsTool(const std::string& type, const std::string& name, const IInterface* parent);

    StatusCode initialize();
    StatusCode calculate();

private: 

    // Methods to fill the Pt variables
    void fillPtFromPtOverlay(Event::PtOverlay*  ptOverlay);
    void recalculatePtVals(Event::EventOverlay* eventOverlay);

    StringProperty          m_root_tree;
    StringArrayProperty     m_pointingHistory;///< history file name and launch date

    astro::PointingHistory* m_history;
    bool                    m_horizontal;

    // Run and event number from overlay event
    unsigned int            m_evtRun;
    unsigned int            m_evtEventId;
    unsigned long long      m_evtEventId64;

    unsigned int            m_triggerBits;

    // ACD hit information
    int                     m_numAcdHits;
    float                   m_acdDepositedEnergy;

    // Cal information
    int                     m_numCalHits;
    float                   m_calDepositedEnergy;

    // Tkr information
    int                     m_tkrStripsHit;

    // Variables to be used to calculate the Pt variables and stored into the ntuple
    double                  m_start;
    float                   m_sc_position[3];
    float                   m_lat_geo;
    float                   m_lon_geo;
    float                   m_lat_mag;
    float                   m_rad_geo;
    float                   m_ra_scz;
    float                   m_dec_scz;
    float                   m_ra_scx; 
    float                   m_dec_scx;
    float                   m_zenith_scz;       ///< space craft zenith angle
    float                   m_B;                ///< magnetic field
    float                   m_L;                ///< McIllwain L parameter

    float                   m_lambda;
    float                   m_R;
    float                   m_bEast;
    float                   m_bNorth;
    float                   m_bUp;
};
//------------------------------------------------------------------------


static const ToolFactory<OverlayValsTool>  Factory;
const IToolFactory& OverlayValsToolFactory = Factory;

//------------------------------------------------------------------------
//! ctor
OverlayValsTool::OverlayValsTool(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent)
                                 : ValBase( type, name, parent ) , m_horizontal(false)
{
    // Declare additional interface
    declareInterface<IValsTool>(this); 

    // declare properties with setProperties calls
    declareProperty("pointing_info_tree_name", m_root_tree = "MeritTuple");

    // doublet, filename and launch date
    declareProperty("PointingHistory",         m_pointingHistory); 
}


/** @page anatup_vars 
    @section Pt  Pt Variables

    These items are added to the merit tuple  to give the current instrument orientation 

<table>
<tr><th> Variable <th>Type<th> Description
<tr><td> OvrEvtRun 
<td>U<td>   Run number, copied from the event header NEW: replaces Run in the merit ntuple  
<tr><td> OvrEvtEventId 
<td>U<td>   Sequence number of event in the run
<tr><td> OvrEvtEventId64
<tr><td> OvrTriggerBits
<td>U<td>   The overlay event's trigger bits
<td>UL<td>   Sequence number of event in the run, 64-bit unsigned version
<tr><td> OvrPtTime       <td>D<td> (s) Current time, same as the elapsed time
<tr><td> OvrPtLat,PtLon  <td>F<td> (deg) lattitude and longitude
<tr><td> OvrPtAlt        <td>F<td> (km) altitude
<tr><td> OvrPtMagLat     <td>F<td> (deg) magnetic latitude, signed; see PtLambda
<tr><td> OvrPtPos[3]     <td>F<td> (m) current orbit position
<tr><td> OvrPtRax,PtDecx <td>F<td> (deg) equatorial coordinates for orientation of S/C x-axis
<tr><td> OvrPtRaz,PtDecz <td>F<td> (deg) equatorial coordinates for orientation of S/C z-axis 
<tr><td> OvrPtSCzenith   <td>F<td> (deg) current angle between zenith and S/C z-axis 
                                   Now signed... + means rocked north
<tr><td> OvrPtMcIlwainB  <td>F<td> McIlwain-L parameter
<tr><td> OvrPtMcIlwainL  <td>F<td> McIlwain-B parameter
<tr><td> OvrPtLambda     <td>F<td> (deg)Lambda parameter, signed according to whether PlLambda 
                                increases (+) or decreases (-) with increasing PtLat
                                Test is done for a one-degree increment of PtLat
<tr><td> OvrPtR          <td>F<td> distance to the dipole center, in Earth radii
<tr><td> OvrPtBEast      <td>F<td> East component of the magnetic field
<tr><td> OvrPtBNorth     <td>F<td> North component of the magnetic field
<tr><td> OvrPtBUp        <td>F<td> Upward component of the magnetic field
</table>
*/

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode OverlayValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the GPS instance: 
    gps = astro::GPS::instance();

    //set the input file to be used as the pointing database, if used
    if( m_pointingHistory.value().empty())
    {
        log << MSG::WARNING << "No history file specified, using default" << endreq;
    }
    else
    {
        std::string filename(m_pointingHistory.value()[0]);
        facilities::Util::expandEnvVar(&filename);
        double offset = 0;
        std::string jStr;
        if( m_pointingHistory.value().size()>1)
        {
            std::string field(m_pointingHistory.value()[1]);
            if(! field.empty() )  // allow null string
            {
                facilities::Timestamp jt(m_pointingHistory.value()[1]);
                offset = (astro::JulianDate(jt.getJulian())
                       - astro::JulianDate::missionStart())
                       * astro::JulianDate::secondsPerDay;
                astro::JulianDate jDate(astro::JulianDate(jt.getJulian()));
                jStr = jDate.getGregorianDate();
            }
        }

        if( m_pointingHistory.value().size() > 2)
        {
            std::string field(m_pointingHistory.value()[2]);
            m_horizontal =! field.empty();
        }

        log << MSG::INFO << "Loading Pointing History File : " << filename <<endreq;

        if( offset>0 )
        {
            log << MSG::INFO   << " with MET offset " ;
            log.precision(12);
            log << offset << " ";
            log.precision(6);
            log << jStr << endreq;
        }

        if( m_horizontal)
        {
            log << MSG::INFO 
                << "   Will override x-direction to be horizontal"<<endreq;
        }

        gps->setPointingHistoryFile(filename, offset, m_horizontal);
    }
   
    // Add the tuple items
    addItem("OvrEvtRun",       &m_evtRun);
    addItem("OvrEvtEventId",   &m_evtEventId);
    addItem("OvrEvtEventId64", &m_evtEventId64);
    addItem("OvrTriggerBits",  &m_triggerBits);
    addItem("OvrNumAcdHits",   &m_numAcdHits);
    addItem("OvrAcdDepEnergy", &m_acdDepositedEnergy);
    addItem("OvrNumCalHits",   &m_numCalHits);
    addItem("OvrCalDepEnergy", &m_calDepositedEnergy);
    addItem("OvrNumTkrHits",   &m_tkrStripsHit);
    addItem("OvrPtTime",       &m_start);
    addItem("OvrPtPos[3]",      m_sc_position);
    addItem("OvrPtLat",        &m_lat_geo);
    addItem("OvrPtLon",        &m_lon_geo);
    addItem("OvrPtMagLat",     &m_lat_mag);
    addItem("OvrPtAlt",        &m_rad_geo);
    addItem("OvrPtRaz",        &m_ra_scz);
    addItem("OvrPtDecz",       &m_dec_scz);
    addItem("OvrPtRax",        &m_ra_scx);
    addItem("OvrPtDecx",       &m_dec_scx);
    addItem("OvrPtSCzenith",   &m_zenith_scz);
    addItem("OvrPtMcIlwainB",  &m_B);
    addItem("OvrPtMcIlwainL",  &m_L);

    addItem("OvrPtLambda",     &m_lambda);
    addItem("OvrPtR",          &m_R);
    addItem("OvrPtBEast",      &m_bEast);
    addItem("OvrPtBNorth",     &m_bNorth);
    addItem("OvrPtBUp",        &m_bUp);

    zeroVals();

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode OverlayValsTool::calculate()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    //
    // Purpose: set tuple items
 
    // Look up the Overlay event timing information
    SmartDataPtr<Event::EventOverlay> header(m_pEventSvc, OverlayEventModel::Overlay::EventOverlay);

    // No overlay event, no variables filled
    if (header)
    {
        // Start with run,event stuff
        m_evtRun         = header->run();
        m_evtEventId     = header->event();
        m_evtEventId64   = header->event();

        // Acd overlay information
        SmartDataPtr<Event::AcdOverlayCol> acdOverlayCol(m_pEventSvc, OverlayEventModel::Overlay::AcdOverlayCol);
        if (acdOverlayCol)
        {
            // Loop through input ACD overlay objects, counting and summing deposited energy
            for(Event::AcdOverlayCol::iterator overIter  = acdOverlayCol->begin(); overIter != acdOverlayCol->end(); overIter++)
            {
                Event::AcdOverlay* acdOverlay = *overIter;

                m_acdDepositedEnergy += acdOverlay->getEnergyDep();
                m_numAcdHits++;
            }
        }

        // Cal Overlay information
        SmartDataPtr<Event::CalOverlayCol> calOverlayCol(m_pEventSvc, OverlayEventModel::Overlay::CalOverlayCol);
        if(calOverlayCol)
        {
            // Loop through input Cal overlay objects, counting and summing deposited energy
            for(Event::CalOverlayCol::iterator overIter  = calOverlayCol->begin(); overIter != calOverlayCol->end(); overIter++)
            {
                Event::CalOverlay* calOverlay = *overIter;

                m_calDepositedEnergy += calOverlay->getEnergy();
                m_numCalHits++;
            }
        }

        // Tkr Overlay information
        SmartDataPtr<Event::TkrOverlayCol> tkrOverlayCol(m_pEventSvc, OverlayEventModel::Overlay::TkrOverlayCol);
        if (tkrOverlayCol)
        {
            // Loop through input Cal overlay objects, counting and summing deposited energy
            for(Event::TkrOverlayCol::iterator overIter  = tkrOverlayCol->begin(); overIter != tkrOverlayCol->end(); overIter++)
            {
                Event::TkrOverlay* tkrOverlay = *overIter;

                // number of hit strips
                m_tkrStripsHit += tkrOverlay->size();
            }
        }

        // Look up the Overlay event timing information
        SmartDataPtr<Event::PtOverlay> ptOverlay(m_pEventSvc, OverlayEventModel::Overlay::PtOverlay);

        if (ptOverlay) fillPtFromPtOverlay(ptOverlay);
        else           recalculatePtVals(header);
 
        // Look up the Overlay event timing information
        SmartDataPtr<Event::GemOverlay> gemOverlay(m_pEventSvc, OverlayEventModel::Overlay::GemOverlay);
        if (gemOverlay) m_triggerBits = gemOverlay->getConditionSummary();
    }
    
    return StatusCode::SUCCESS;
}

void OverlayValsTool::fillPtFromPtOverlay(Event::PtOverlay*  ptOverlay)
{
    m_start          = ptOverlay->getStartTime();
    m_sc_position[0] = ptOverlay->getSC_Position()[0];
    m_sc_position[1] = ptOverlay->getSC_Position()[1];
    m_sc_position[2] = ptOverlay->getSC_Position()[2];
    m_lat_geo        = ptOverlay->getLatGeo();
    m_lon_geo        = ptOverlay->getLonGeo();
    m_lat_mag        = ptOverlay->getLatMag();
    m_rad_geo        = ptOverlay->getRadGeo();
    m_ra_scz         = ptOverlay->getRaScz();
    m_dec_scz        = ptOverlay->getDecScz();
    m_ra_scx         = ptOverlay->getRaScx(); 
    m_dec_scx        = ptOverlay->getDecScx();
    m_zenith_scz     = ptOverlay->getZenithScz(); 
    m_B              = ptOverlay->getB();
    m_L              = ptOverlay->getL();
    m_lambda         = ptOverlay->getLambda();
    m_R              = ptOverlay->getR();
    m_bEast          = ptOverlay->getBEast();
    m_bNorth         = ptOverlay->getBNorth();
    m_bUp            = ptOverlay->getBUp();

    return;
}

void OverlayValsTool::recalculatePtVals(Event::EventOverlay* eventOverlay)
{
    // get event time from header or merit 
    // and look up position info from the history  

    double etime = eventOverlay->time();

    gps->time(etime);

    m_start = gps->time();
    CLHEP::Hep3Vector pos_km = gps->position();
    CLHEP::Hep3Vector location = 1.e3* pos_km; 

    // cartesian location of the LAT (in m)
    m_sc_position[0] = location.x();
    m_sc_position[1] = location.y();
    m_sc_position[2] = location.z(); 

    m_ra_scx =    gps->xAxisDir().ra();
    m_dec_scx =   gps->xAxisDir().dec();

    m_ra_scz =    gps->zAxisDir().ra();
    m_dec_scz =   gps->zAxisDir().dec();

    const astro::EarthCoordinate& earth(gps->earthpos());

    m_lat_geo = earth.latitude(); 
    m_lon_geo = earth.longitude(); 

    m_rad_geo = earth.altitude(); 
    m_L       = earth.L();
    m_B       = earth.B();

    m_R       = earth.R();

    m_lambda  = earth.lambda();

    CLHEP::Hep3Vector bField = earth.magnetic_field();
    m_bEast  = bField.x();
    m_bNorth = bField.y();
    m_bUp    = bField.z();

    m_lat_mag = earth.geolat();
    //m_in_saa= earth.insideSAA()? 1:0;
    m_zenith_scz = gps->zenithDir().difference(gps->zAxisDir());

    // sign the rocking angle
    double temp = pos_km[2]/pos_km.mag();
    if(m_dec_scz < temp) m_zenith_scz *= -1.0;
    m_zenith_scz *= R2D; // and convert to degrees

    return;
}
