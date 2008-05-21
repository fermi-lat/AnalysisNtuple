/** @file ObfCoordsAlg.cxx
@brief Declaration and implementation of Gaudi algorithm ObfCoordsAlg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/ObfCoordsAlg.cxx,v 1.11.4.1 2008/04/28 21:57:21 usher Exp $
*/
// Include files

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "astro/GPS.h"
#include "geometry/Vector.h"
#include "FluxSvc/IFluxSvc.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <cassert>
#include <map>

// forward declatation of the worker
class ObfCworker;

namespace { // anonymous namespace for file-global
    std::string treename("MeritTuple");
    astro::GPS* gps(0);
#include "Item.h"
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class ObfCoordsAlg
@brief Extract info from tuple, etc. to add ObfCoords items to this of another tree
*/
class ObfCoordsAlg : public Algorithm {

public:
    ObfCoordsAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
private:
    /// this guy does the work!
    ObfCworker * m_worker;
    //counter
    int m_count;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static const AlgFactory<ObfCoordsAlg>  Factory;
const IAlgFactory& ObfCoordsAlgFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ObfCworker{ 
public:
    ObfCworker();

    void evaluate();

private:
    // tuple items expect to find
    //TypedItem<unsigned int, 'i'> EvtRun, EvtEventId;

    // These will replace Filter%Dir's above
    Item GrbXDir;
    Item GrbYDir;
    Item GrbZDir;

    //ObfCoords entries to create
    float m_obfRa,m_obfDec,m_obfL,m_obfB;
    float m_grbRa,m_grbDec,m_grbL,m_grbB;
    //float m_obfZen,m_obfAzim;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ObfCoordsAlg::ObfCoordsAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator)
, m_count(0)
{
    declareProperty("TreeName",  treename="MeritTuple");
}

StatusCode ObfCoordsAlg::initialize()
{
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // Use the Job options service to get the Algorithm's parameters
    setProperties();

    // get a pointer to RootTupleSvc 
    if( (sc = service("RootTupleSvc", rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
        return sc;
    }
    m_worker = new ObfCworker();

    // get the GPS instance: either from FluxSvc or local, non-MC mode
    IFluxSvc* fluxSvc(0);
    if( service("FluxSvc", fluxSvc).isFailure() ){
        //log << MSG::INFO << "Will use the GPS singleton to get coordinates" << endreq;
        gps = astro::GPS::instance();
    }else{
        gps = fluxSvc->GPSinstance();
    }
    return sc;
}

StatusCode ObfCoordsAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    m_count++;

    // now have the worker do it
    m_worker->evaluate();

    return sc;
}

StatusCode ObfCoordsAlg::finalize()
{
    return StatusCode::SUCCESS;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/** @page anatup_vars 

@section ObfCoords  Celestial Coordinates of the Best OnboardFilter Track

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> FilterRa, FilterDec 
<td>F<td>  (deg) reconstructed direction in equatorial coordinates       
<tr><td> FilterL, FilterB 
<td>F<td>  (deg) galactic longitude and latitude of reconstructed direction
</table> 
*/


ObfCworker::ObfCworker()
// initialize pointers to current items
: GrbXDir("GrbXDir")
, GrbYDir("GrbYDir")
, GrbZDir("GrbZDir")
{
    //now create new items 

    addItem( "FilterRa",  m_obfRa);
    addItem( "FilterDec", m_obfDec);
    addItem( "FilterL",   m_obfL);
    addItem( "FilterB",   m_obfB);

    addItem( "GrbRa",     m_grbRa);
    addItem( "GrbDec",    m_grbDec);
    addItem( "GrbL",      m_grbL);
    addItem( "GrbB",      m_grbB);

    //addItem( "ObfcZenithTheta",   m_obfZen);
    //addItem( "ObfcEarthAzimuth",  m_obfAzim);
}


void ObfCworker::evaluate()
{

    m_obfRa = m_obfDec = m_obfL = m_obfB = 0;
    m_grbRa = m_grbDec = m_grbL = m_grbB = 0;
    // convert to (ra, dec)

    // "Best" GRB track in the Grb tuple variables
    Vector grbDir(GrbXDir, GrbYDir, GrbZDir);
    if (grbDir.mag()==0) return;
    astro::SkyDir skyGrbdir( gps->toSky(-grbDir) );
    m_grbRa  = skyGrbdir.ra();
    m_grbDec = skyGrbdir.dec();
    m_grbL   = skyGrbdir.l();
    m_grbB   = skyGrbdir.b();

    // Old school stuff first
    Vector filterDir(GrbXDir, GrbYDir, GrbZDir);
    if (filterDir.mag()==0) return;
    // Filter direction points up... 
    // toSky converts a *particle* direction
    // into a direction on the sky, so the minus-sign is needed below (twice)!
    astro::SkyDir skydir( gps->toSky(-filterDir) );
    m_obfRa  = skydir.ra();
    m_obfDec = skydir.dec();
    m_obfL   = skydir.l();
    m_obfB   = skydir.b();

    return;
}
//------------------------------------------------------------------------------

