/** @file GltValsTool.cxx
@brief Calculates the Trigger analysis variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/GltValsTool.cxx,v 1.13 2004/09/30 17:12:21 lsrea Exp $
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
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

/*! @class GltValsTool
@brief calculates trigger values

@authors Bill Atwood, Leon Rochester
*/

namespace {
    const int _nTowers = 16; // maximum possible number of towers
    const int _nLayers = 18; // maximum possible number of layers

    //Replace these with official enum when it exists
    enum { 
        ACDLO = 0, 
        ACDHI = 1, 
        TRACK = 2, 
        CALLO = 3, 
        CALHI = 4, 
        THROTTLE = 5
    }; 
    enum { SHIFTGem = 8 };
    enum { MASKGlt = 0xff, MASKGem = 0xff } ;
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

private:
    /// pointer to tracker geometry
    ITkrGeometrySvc*       pTkrGeoSvc;

    //TkrClusters Tuple Items
    double Trig_word;
    double Trig_GemSummary; // new
    double Trig_tower;
    double Trig_xTower;
    double Trig_yTower; 
    double Trig_layer;
    double Trig_total; 
    double Trig_numTowers;
    double Trig_type; // was: 1= corner, 2 = side, 3 = core, now number of exposed sides (0-4)
    double Trig_moment; 
    double Trig_zDir; 

    ITkrQueryClustersTool* m_clusTool;
};

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
        if(service( "TkrGeometrySvc", pTkrGeoSvc, true ).isFailure()) {
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
    addItem("GltGemSummary", &Trig_GemSummary); // new
    addItem("GltTower",      &Trig_tower); 
    addItem("GltXTower",     &Trig_xTower);
    addItem("GltYTower",     &Trig_yTower);
    addItem("GltLayer",      &Trig_layer); 
    addItem("GltTotal",      &Trig_total);
    addItem("GltNumTowers",  &Trig_numTowers);
    addItem("GltType",       &Trig_type);  
    addItem("GltMoment",     &Trig_moment);
    addItem("GltZDir",       &Trig_zDir);  

    zeroVals();

    return sc;
}


StatusCode GltValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    // m_pEventSvc alreay checked by doCalcIfNotDone, no need to repeat

    // Recover EventHeader Pointer
    SmartDataPtr<Event::EventHeader> 
        pEvent(m_pEventSvc, EventModel::EventHeader);
    // Recover Track associated info.

    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    //not currently used:
    //SmartDataPtr<Event::TkrFitTrackCol>    
    //    pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);

    int iTrig_tower  = -1;
    int iTrig_layer  = 18;
    int xTower = -1;
    int yTower = -1; 
    int iTrig_type = 0;

    int nLayers  = pTkrGeoSvc->numLayers();
    int nXTowers = pTkrGeoSvc->numXTowers();
    int nYTowers = pTkrGeoSvc->numYTowers();
    int nTowers  = nXTowers*nYTowers;

    if(!pEvent || !pClusters) return StatusCode::FAILURE;

    unsigned int word = pEvent->trigger();
    // Waiting for official Trigger enum list

    // This is the same as the old GltWord
    Trig_word = word & MASKGlt;
    Trig_GemSummary = (word >> SHIFTGem) & MASKGem;

    bool three_in_a_row = ((word & (1<<TRACK))!=0);

    int tower, layer;
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
            Event::TkrClusterVec xHitList = m_clusTool->getClustersReverseLayer(idents::TkrId::eMeasureX,layer);
            Event::TkrClusterVec yHitList = m_clusTool->getClustersReverseLayer(idents::TkrId::eMeasureY,layer);

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
        // If more than one tower triggers, use the one with the lowest recon layer number
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
                    if(layer < iTrig_layer) {
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
            iTrig_type = pTkrGeoSvc->getTowerType(iTrig_tower);
        }

        // Now find the average location of all hits
        // This is probably too coarse to be useful
        double x_sum, y_sum, z_sum, wts, g1, g2, L11, L12, L22, zmax, zmin;
        x_sum = y_sum = z_sum = wts = g1 = g2 = L11 = L12 = L22 = 0.;
        zmax = 0.;
        zmin = 20.; 
        HepSymMatrix moments(3,0); 

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
                HepMatrix U(3,3);
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
    Trig_xTower = xTower;
    Trig_yTower = yTower;
    Trig_type   = iTrig_type;

    return sc;
}
