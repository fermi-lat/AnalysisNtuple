/** @file GltValsTool.cxx
@brief Calculates the Trigger analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/GltValsTool.cxx,v 1.41 2008/09/10 01:24:33 lsrea Exp $
*/

// Include files

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/GaudiException.h" 

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "idents/TowerId.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "LdfEvent/EventSummaryData.h"
#include "LdfEvent/LsfMetaEvent.h"
#include "LdfEvent/Gem.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "enums/TriggerBits.h"

/*! @class GltValsTool
@brief calculates trigger values

@authors Bill Atwood, Leon Rochester
*/

namespace {
    const int _nTowers = 16; // maximum possible number of towers
    const int _nLayers = 18; // maximum possible number of layers
    const unsigned int _unset   = (unsigned int) -1; // default for some of the variables
}

class GltValsTool : public ValBase
{
public:

    GltValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~GltValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

    void zeroVals();

private:
    /// pointer to tracker geometry
    ITkrGeometrySvc*       m_tkrGeom;

    //TkrClusters Tuple Items
    float Trig_word;
    float Trig_GemSummary;
    float Trig_evtFlags;
    float Trig_tower;
    float Trig_xTower;
    float Trig_yTower; 
    float Trig_layer;
    float Trig_total; 
    float Trig_numTowers;
    float Trig_type; // was: 1= corner, 2 = side, 3 = core, now number of exposed sides (0-4)
    float Trig_moment; 
    float Trig_zDir; 
    int   Trig_engine;
    int   Trig_gemengine;
    int   Trig_gltprescale;
    int   Trig_gemprescale;
    int   Trig_prescaleexpired;
    int   Trig_sourcegps;
    int   Trig_gemDeltaEventTime;
    int   Trig_gemDeltaWindowOpenTime;
    int   Trig_eventSize;
    int   Trig_compressedEventSize;
    unsigned long long Trig_elapsedTime;
    unsigned long long Trig_liveTime;

    ITkrQueryClustersTool* m_clusTool;
};

/** @page anatup_vars 
@section gltvalstool GltValsTool Variables
(Glt = GLAST Trigger)

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> GltWord 
<td>F<td>   Trigger word from the event header, based on digi information (null value: -1)
The definition of the bits can be found in TriggerBits.h in the enums package.

@verbatim
b_ROI  =    0x01  ///>  Trottle bit
b_Track=    0x02  ///>  3 consecutive x-y layers hit
b_LO_CAL=   0x04  ///>  single log above low threshold
b_HI_CAL=   0x08  ///>  single log above high threshold
b_CNO  =    0x10  ///>  cover or side veto, high threshold ("CNO")
b_trottle=  0x20  ///>  temporary copy of the ROI bit, for backward compatibilty
b_ACDL =    0x40  ///> set if cover or side veto, low threshold
@endverbatim
<tr><td> GltGemSummary 
<td>F<td>   For Data, a copy of the GEM condtion summary; for MC, 
the same as calculated from the trigger (digi) information (null value: -1) 
The bits are defined in GemConditionsSummary.h in the enums package.

Here are the definitions as of GlastRelease v7r2 (They are not expected to change.):
@verbatim
ROI   = 0x01         CNO       =  0x10
TKR   = 0x02         PERIODIC  =  0x20
CALLE = 0x04         SOLICITED =  0x40
CALHE = 0x08         EXTERNAL  =  0x80
@endverbatim
<tr><td>  GltEventFlags  
<td>F<td>   Event flags from Ebf 
<tr><td> GltTower 
<td>F<td>   Id of tower that triggered; if more than one, 
choose the tower with the highest triggered layer (null value: -1) 
<tr><td> Glt[X/Y]Tower 
<td>F<td>   [x/y] index of GltTower (null value = -1) 
<tr><td> GltTotal 
<td>F<td>   Number of possible triggers (a straight-through track in a single tower 
produces 13 potential triggers 
<tr><td> GltNumTowers 
<td>F<td>   Number of towers which trigger (crossing track) 
<tr><td> GltType 
<td>F<td>   Number of exposed sides for the triggered tower: <br>
            0 = central tower, 1 = side tower, <br> 
            2 = edge edge tower, 4 = any single-tower setup 	
<tr><td> GltMoment 
<td>F<td>   Do not use 
<tr><td> GltZDir 
<td>F<td>   Do not use 
<tr><td> GltEngine 
<td>I<td>   The engine number corresponding to GltWord  
<tr><td> GltGemEngine
<td>I<td>   The engine number corresponding to the GltGemSummary
<tr><td> GltEnginePrescale
<td>I<td>   The GLT prescale factor for the event
<tr><td> GltGemEnginePrescale
<td>I<td>   The GEM prescale factor for the event
<tr><td> GltPrescaleExpired
<td>I<td>   true if the prescale counter expired for this event
<tr><td> GltSourceGps
<td>I<td>   true if timing comes from the GPS
<tr><td> GltGemDeltaEventTime
<td>I<td>   GEM Delta event time counter in ticks of 50 ns. 
            Saturates at 3.3 ms
<tr><td> GltGemDeltaWindowOpenTime
<td>I<td>   GEM Delta window open time counter: Number of
            of 50 ns ticks between the opening of the trigger window 
            of the previous event and of this event. 
            Wraps around at ~3.3 msec
<tr><td> GltEventSize
<td>I<td>   Uncompressed size of the event in bytes
<tr><td> GltCompressedEventSize
<td>I<td>   Compressed size of the event in bytes
<tr><td> GltGemElapsedTime
<td>UL<td>   Gem raw count scaler (64-bit, counts at 20 MHz).
<tr><td> GltGemLiveTime
<td>UL<td>   Gem livetime scaler (64-bit, counts at 20 MHz),
             increments every 50-nsec tick when instrument is live
            </table>
*/


// Static factory for instantiation of algtool objects
static ToolFactory<GltValsTool> s_factory;
const IToolFactory& GltValsToolFactory = s_factory;

// Standard Constructor
GltValsTool::GltValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

StatusCode GltValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the services
    StatusCode fail = StatusCode::FAILURE;
    if( serviceLocator() ) {       
        if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return fail;
        }
    } else {
        return fail;
    }

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }

    // load up the map

    addItem("GltWord",       &Trig_word);
    addItem("GltGemSummary", &Trig_GemSummary);
    addItem("GltEventFlags", &Trig_evtFlags);    // new
    addItem("GltTower",      &Trig_tower); 
    addItem("GltXTower",     &Trig_xTower);
    addItem("GltYTower",     &Trig_yTower);
    addItem("GltLayer",      &Trig_layer); 
    addItem("GltTotal",      &Trig_total);
    addItem("GltNumTowers",  &Trig_numTowers);
    addItem("GltType",       &Trig_type);  
    addItem("GltMoment",     &Trig_moment);
    addItem("GltZDir",       &Trig_zDir);  
    addItem("GltEngine",     &Trig_engine);  
    addItem("GltGemEngine",  &Trig_gemengine);
    addItem("GltGemDeltaEventTime",      &Trig_gemDeltaEventTime);
    addItem("GltGemDeltaWindowOpenTime", &Trig_gemDeltaWindowOpenTime);
    addItem("GltEnginePrescale",         &Trig_gltprescale);
    addItem("GltGemEnginePrescale",      &Trig_gemprescale);
    addItem("GltPrescaleExpired",        &Trig_prescaleexpired);
    addItem("GltSourceGps",              &Trig_sourcegps);

    addItem("GltEventSize",              &Trig_eventSize);
    addItem("GltCompressedEventSize",    &Trig_compressedEventSize);

    addItem("GltGemElapsedTime",         &Trig_elapsedTime);
    addItem("GltGemLiveTime",            &Trig_liveTime);

    zeroVals();

    return sc;
}

void GltValsTool::zeroVals()
{
    ValBase::zeroVals();

    Trig_gemprescale     = _unset;
    Trig_gltprescale     = _unset;
    Trig_prescaleexpired = _unset;
    Trig_gemDeltaEventTime      = -1;
    Trig_gemDeltaWindowOpenTime = -1;
}

StatusCode GltValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    // m_pEventSvc already checked by doCalcIfNotDone, no need to repeat

    int nLayers  = m_tkrGeom->numLayers();
    int nXTowers = m_tkrGeom->numXTowers();
    int nYTowers = m_tkrGeom->numYTowers();
    int nTowers  = nXTowers*nYTowers;


    int iTrig_tower  = -1;
    int iTrig_layer  = -1;
    int iTrig_xTower = -1;
    int iTrig_yTower = -1; 
    int iTrig_type   =  0;

    Trig_layer  = iTrig_layer;
    Trig_tower  = iTrig_tower;
    Trig_xTower = iTrig_xTower;
    Trig_yTower = iTrig_yTower;
    Trig_type   = iTrig_type;

    // Recover EventHeader Pointer
    SmartDataPtr<Event::EventHeader> 
        pEvent(m_pEventSvc, EventModel::EventHeader);
    unsigned int word  = 0;
    unsigned int word2 = _unset;
    
    if(pEvent) {
    word  = pEvent->trigger();
    word2 = pEvent->triggerWordTwo();
     // by default, engines are initialized to "_unset"; see zeroVals() above
    Trig_gemprescale     = pEvent->gemPrescale();
    Trig_gltprescale     = pEvent->gltPrescale();
    Trig_prescaleexpired = pEvent->prescaleExpired();
    }

    // This is the same as the old GltWord
    // actually only 6 bits, but no harm (I think!)  
    // note that the trigger word from the header has 3 GEM_offset fields that we unpack here:
    Trig_word = word & enums::GEM_mask;  // the GltWord, set for simulation from hits or digis

    Trig_GemSummary = (word >> enums::GEM_offset) & enums::GEM_mask; // the GEM condition word, same as previous if simulation
    
    Trig_engine = (word >> (2*enums::GEM_offset)) & enums::GEM_mask; // the trigger engine number

    unsigned int Trig_gltengine = word2 & enums::ENGINE_mask;
    // If the engine number is set, then we will store this one from the TrgConfigSvc in the tuple
    if (Trig_gltengine != enums::ENGINE_unset) Trig_engine = Trig_gltengine;
    Trig_gemengine = ((word2 >> enums::ENGINE_offset) & enums::ENGINE_mask); // GEM trigger engine number

    SmartDataPtr<LsfEvent::MetaEvent> 
        metaTds(m_pEventSvc, "/Event/MetaEvent");
    if(metaTds) {
        Trig_sourcegps =  metaTds->time().current().sourceGps();
        Trig_compressedEventSize = metaTds->compressedSize();
        Trig_elapsedTime = metaTds->scalers().elapsed();
        Trig_liveTime    = metaTds->scalers().livetime();
    }

    SmartDataPtr<LdfEvent::EventSummaryData> 
        eventSummary(m_pEventSvc, "/Event/EventSummary"); 
    if(eventSummary) {
        Trig_evtFlags =  eventSummary->eventFlags();
        Trig_eventSize = eventSummary->eventSizeInBytes();
    }

    // GemDeltaTimes
    SmartDataPtr<LdfEvent::Gem> gemTds(m_pEventSvc, "/Event/Gem");
    if(gemTds) {
        Trig_gemDeltaEventTime      = gemTds->deltaEventTime();
        Trig_gemDeltaWindowOpenTime = gemTds->deltaWindowOpenTime();
    }

    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    // everything from here on out uses Clusters
    if(!pClusters) return sc;

    bool three_in_a_row = ((word & enums::b_Track)!=0);

    int tower, layer;
    // needs to be recast into an indexed vector maybe
    short x_hits[_nTowers][_nLayers]; 
    short y_hits[_nTowers][_nLayers]; 
    for(tower=0; tower<_nTowers; ++tower) {
        for(layer=0; layer<_nLayers; ++layer){
            x_hits[tower][layer]=0;
            y_hits[tower][layer]=0;
        }
    }

    //Search for x-y paired layer.... 
    if (pClusters && three_in_a_row)
    {
        // Make a hit count per plane, per tower.  
        layer = nLayers;
        while(layer--)
        {
            Event::TkrClusterVec xHitList = 
                m_clusTool->getClusters(idents::TkrId::eMeasureX,layer);
            Event::TkrClusterVec yHitList = 
                m_clusTool->getClusters(idents::TkrId::eMeasureY,layer);

            int x_hitCount = xHitList.size(); 
            int y_hitCount = yHitList.size();
            if(x_hitCount > 0 && y_hitCount > 0) {

                int hit;
                // x hits
                for(hit=0; hit<x_hitCount; ++hit){
                    tower = xHitList[hit]->tower();
                    x_hits[tower][layer]++; 
                }
                // y hits
                for(hit=0; hit<y_hitCount; ++hit){
                    tower = yHitList[hit]->tower();
                    y_hits[tower][layer]++;
                }
            }
        }

        // Now search for the 3-in-a-rows (need to look for x & y hits in each...)
        // If more than one tower triggers, use the one with the highest layer number
        int nTotTrigs = 0;   // total number of possible triggers
        int nTowerTrigs = 0; // total number of towers that could trigger
        for(tower = 0; tower<nTowers; ++tower) {
            int nTrigsThis = 0;
            for(layer=0; layer<nLayers-2; ++layer){
                if( (x_hits[tower][layer]   > 0 && y_hits[tower][layer]   > 0 ) &&
                    (x_hits[tower][layer+1] > 0 && y_hits[tower][layer+1] > 0 ) &&
                    (x_hits[tower][layer+2] > 0 && y_hits[tower][layer+2] > 0 )) 
                {
                    nTrigsThis++;
                    if(layer > iTrig_layer) {
                        iTrig_layer = layer;
                        iTrig_tower = tower;
                    }
                }
            }
            nTotTrigs += nTrigsThis;
            if(nTrigsThis>0) nTowerTrigs++;
        }
        Trig_total  = nTotTrigs;
        Trig_numTowers = nTowerTrigs;

        // Now classify according to tower type
        // new classification is number of exposed sides
        if(iTrig_tower >= 0) {
            iTrig_type = m_tkrGeom->getTowerType(iTrig_tower);
            idents::TowerId towerId = idents::TowerId(iTrig_tower);
            iTrig_xTower = towerId.ix();
            iTrig_yTower = towerId.iy();
        }

        // Now find the average location of all hits
        // This is probably too coarse to be useful
        double x_sum, y_sum, z_sum, wts, g1, g2, L11, L12, L22, zmax, zmin;
        x_sum = y_sum = z_sum = wts = g1 = g2 = L11 = L12 = L22 = 0.;
        zmax = 0.;
        zmin = 20.; 
        CLHEP::HepSymMatrix moments(3,0); 

        if(iTrig_tower >= 0) {
            for(tower = 0; tower<_nTowers; ++tower) {
                idents::TowerId towerId = idents::TowerId(tower);
                int ix = towerId.ix();
                int iy = towerId.iy();
                if (ix+1==nXTowers || iy+1==nYTowers) continue;
                for(layer=0; layer<nLayers; ++layer){
                    int x_counts = x_hits[tower][layer];
                    int y_counts = y_hits[tower][layer];
                    double weight = x_counts+y_counts;
                    if(x_counts>0 && y_counts>0) { // Don't count random noise hits
                        wts += weight;
                        double x = (ix+1)*5.;
                        double y = (iy+1)*5.;
                        double z = (layer+1);
                        x_sum += x * weight;
                        y_sum += y * weight;
                        z_sum += z * weight;
                        g1    += 1.;
                        g2    += z;
                        L11   += 1./weight;
                        L12   += z /weight;
                        L22   += z*z/weight; 
                        if(zmax < z ) zmax = z;
                        if(zmin > z ) zmin = z; 

                    }
                }
            }

            // Now form moment tensor, etc..... 
            if(wts > 0) {
                x_sum /= wts;
                y_sum /= wts;
                z_sum /= wts;
                for(tower = 0; tower<_nTowers; ++tower) {
                    idents::TowerId towerId = idents::TowerId(tower);
                    int ix = towerId.ix();
                    int iy = towerId.iy();
                    if (ix+1==nXTowers || iy+1==nYTowers) continue;
                    for(layer=0; layer<nLayers-2; ++layer){
                        int x_counts = x_hits[tower][layer];
                        int y_counts = y_hits[tower][layer];
                        double weight = x_counts+y_counts;
                        if(x_counts>0 && y_counts>0) { // Don't count random noise hits
                            double xres = (ix+1)*5.  - x_sum;
                            double yres = (iy+1)*5.  - y_sum;
                            double zres = (layer+1)  - z_sum;
                            moments(1,1) += xres * xres * weight;
                            moments(1,2) += xres * yres * weight;
                            moments(1,3) += xres * zres * weight;
                            moments(2,2) += yres * yres * weight;
                            moments(2,3) += yres * zres * weight;
                            moments(3,3) += zres * zres * weight;
                        }
                    }
                }
                moments /= wts;

                // Now find transformation to diagonalize
                CLHEP::HepMatrix U(3,3);
                U = diagonalize(&moments);

                // Now get directions, etc. 
                double Det = L11*L22 - L12*L12;
                if(fabs(Det) > 0) {
                    double z_hit_slope = (g2*L11 - g1*L12)/Det;
                    Trig_zDir = z_hit_slope*(zmax-zmin)*wts; 
                }
                int sm = 1;
                if(moments(2,2) > moments(1,1))   sm = 2;
                if(moments(3,3) > moments(sm,sm)) sm = 3;

                Trig_moment = moments(sm,sm);
                //            Vector t_moment(U(1,sm),U(2,sm),U(3,sm));
            }
        }
    }

    Trig_layer  = iTrig_layer;
    Trig_tower  = iTrig_tower;
    Trig_xTower = iTrig_xTower;
    Trig_yTower = iTrig_yTower;
    Trig_type   = iTrig_type;

    return sc;
}
