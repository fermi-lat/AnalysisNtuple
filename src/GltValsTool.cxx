// $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/GltValsTool.cxx,v 1.3 2003/02/26 17:38:48 lsrea Exp $

// Include files


#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

namespace {
    const int _nLayers = 18;
    const int _nTowers = 16;
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
    
    //TkrClusters Tuple Items
    double Trig_word;
    double Trig_tower;
    double Trig_xTower;
    double Trig_yTower; 
    double Trig_layer;
    double Trig_total; 
    double Trig_type; // 1= corner, 2 = side, 3 = core
    double Trig_moment; 
    double Trig_zDir; 
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
    
    if( serviceLocator() ) {      
    } else {
        return StatusCode::FAILURE;
    }

    // load up the map

    addItem("TRG_Word",      &Trig_word);  
    addItem("TRG_Tower",     &Trig_tower); 
    addItem("TRG_XTower",    &Trig_xTower);
    addItem("TRG_YTower",    &Trig_yTower);
    addItem("TRG_Layer",     &Trig_layer); 
    addItem("TRG_Total",     &Trig_total); 
    addItem("TRG_Type",      &Trig_type);  
    addItem("TRG_Moment",    &Trig_moment);
    addItem("TRG_zDir",      &Trig_zDir);  
    
    zeroVals();
    
    return sc;
}


StatusCode GltValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    // Recover EventHeader Pointer
    SmartDataPtr<Event::EventHeader> pEvent(m_pEventSvc, EventModel::EventHeader);
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>    
        pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    
    unsigned int word = pEvent->trigger();
    if(word > 1024) return sc; 
    bool three_in_a_row = (word/4)%2 == 1; 
    Trig_word = word;
    
    int tower, layer;
    
    int lat_hits[_nTowers][_nLayers]; 
    for(tower=0; tower<_nTowers; tower++) {
        for(layer=0; layer<_nLayers; layer++){
            lat_hits[tower][layer]=0;
        }
    }
    
    //Search for x-y paired layer.... 
    if (pClusters && three_in_a_row)
    {
        // Make a hit count per plane, per tower.  
        //    X hits 0 - 999 and Y hits 1000 -...
        layer = _nLayers;
        while(layer--)
        {
            int x_hitCount = pClusters->nHits(Event::TkrCluster::X, layer); 
            int y_hitCount = pClusters->nHits(Event::TkrCluster::Y, layer);
            if(x_hitCount > 0 && y_hitCount > 0) {
                std::vector<Event::TkrCluster*> xHitList;
                std::vector<Event::TkrCluster*> yHitList;
                xHitList = pClusters->getHits(Event::TkrCluster::X,layer);
                yHitList = pClusters->getHits(Event::TkrCluster::Y,layer);
                
                for(int ix=0; ix<x_hitCount; ix++){
                    int x_tower = xHitList[ix]->tower();
                    lat_hits[x_tower][layer] += 1; 
                }
                
                for(int iy=0; iy<y_hitCount; iy++){
                    int y_tower = yHitList[iy]->tower();
                    lat_hits[y_tower][layer] += 1000;
                }
            }
        }
        
        // Now search for the 3-in-a-rows (need to look for x & y hits in each...)
        for(tower = 0; tower<_nTowers; tower++) {
            for(layer=0; layer<_nLayers-2; layer++){
                if((lat_hits[tower][layer]%1000   > 0 && lat_hits[tower][layer]/1000   > 0 ) &&
                    (lat_hits[tower][layer+1]%1000 > 0 && lat_hits[tower][layer+1]/1000 > 0 ) &&
                    (lat_hits[tower][layer+2]%1000 > 0 && lat_hits[tower][layer+2]/1000 > 0 )) {
                    Trig_total += 1;
                    if(Trig_layer > layer) {
                        Trig_layer = layer;
                        Trig_tower = tower;
                    }
                }
            }
        }
        
        // Now classify according to tower type
        if(Trig_tower >= 0) {
            int towerId = (int) Trig_tower;
            Trig_xTower = towerId%4;
            Trig_yTower = towerId/4; 
            if(Trig_tower==0  || Trig_tower==3 || 
                Trig_tower==12 || Trig_tower==15) Trig_type = 1;
            else if(Trig_tower==5 || Trig_tower==6 || 
                Trig_tower==9 || Trig_tower==10) Trig_type = 3;
            else Trig_type = 2;
        }
        
        // Now find the average location of all hits
        double x_sum, y_sum, z_sum, wts, g1, g2, L11, L12, L22, zmax, zmin;
        x_sum = y_sum = z_sum = wts = g1 = g2 = L11 = L12 = L22 = 0.;
        zmax = 0.;
        zmin = 20.; 
        HepSymMatrix moments(3,0); 
        
        if(Trig_tower >= 0) {
            for(tower = 0; tower<_nTowers; tower++) {
                int ix = tower%4;
                int iy = tower/4;
                for(layer=0; layer<_nLayers; layer++){
                    int x_counts = lat_hits[tower][layer]%1000;
                    int y_counts = lat_hits[tower][layer]/1000;
                    double weight = x_counts+y_counts;
                    if(weight > 1) { // Don't count random noise hits
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
                for(tower = 0; tower<_nTowers; tower++) {
                    int ix = tower%4;
                    int iy = tower/4;
                    for(layer=0; layer<_nLayers-2; layer++){
                        int x_counts = lat_hits[tower][layer]%1000;
                        int y_counts = lat_hits[tower][layer]/1000;
                        double weight = x_counts+y_counts;
                        if(weight > 1) { // Don't count random noise hits
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
    
    return sc;
}
