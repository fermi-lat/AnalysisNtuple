/** @file McCoordsAlg.cxx
@brief Declaration and implementation of Gaudi algorithm McCoordsAlg

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/McCoordsAlg.cxx,v 1.3 2007/03/30 12:53:48 lsrea Exp $
*/
// Include files

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "astro/GPS.h"
#include "geometry/Vector.h"
#include "FluxSvc/IFluxSvc.h"

//#include "ntupleWriterSvc/INTupleWriterSvc.h"


#include <cassert>
#include <map>

// forward declatation of the worker
class McCworker;

namespace { // anonymous namespace for file-global
    IFluxSvc* fluxSvc;
    std::string treename("MeritTuple");
#include "Item.h"
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class McCoordsAlg
@brief Extract info from tuple, etc. to add McCoords items to this of another tree
*/
class McCoordsAlg : public Algorithm {

public:
    McCoordsAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
private:
    /// this guy does the work!
    McCworker * m_worker;
    //counter
    int m_count;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static const AlgFactory<McCoordsAlg>  Factory;
const IAlgFactory& McCoordsAlgFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class McCworker{ 
public:
    McCworker();

    void evaluate();

private:
    // tuple items expect to find
    //TypedItem<unsigned int, 'i'> EvtRun, EvtEventId;

    // these all float or double
    Item McXDir;
    Item McYDir;
    Item McZDir;

    //McCoords entries to create
    float m_mcRa,m_mcDec,m_mcL,m_mcB;
    float m_mcZen,m_mcAzim;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
McCoordsAlg::McCoordsAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)
{
    declareProperty("TreeName",  treename="MeritTuple");
}

StatusCode McCoordsAlg::initialize()
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
    m_worker = new McCworker();

    service("FluxSvc", fluxSvc, true);

    m_count = 0;

    return sc;
}

StatusCode McCoordsAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    m_count++;

    // now have the worker do it
    m_worker->evaluate();

    return sc;
}

StatusCode McCoordsAlg::finalize()
{
    return StatusCode::SUCCESS;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/** @page anatup_vars 

@section McCoords  Celestial Coordinates of the Primary MC Particle

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> McRa, McDec 
<td>F<td>  (deg) reconstructed direction in equatorial coordinates       
<tr><td> McZenithTheta, McEarthAzimuth 
<td>F<td>  (deg) reconstucted direction with respect to local zenith system
<tr><td> McL, McB 
<td>F<td>  (deg) galactic longitude and latitude of reconstructed direction
</table> 
*/

McCworker::McCworker()
// initialize pointers to current items
: McXDir("McXDir")
, McYDir("McYDir")
, McZDir("McZDir")
{
    //now create new items 

    addItem( "McRa",            m_mcRa);
    addItem( "McDec",           m_mcDec);
    addItem( "McL",             m_mcL);
    addItem( "McB",             m_mcB);
    addItem( "McZenithTheta",   m_mcZen);
    addItem( "McEarthAzimuth",  m_mcAzim);
}

void McCworker::evaluate()
{
    // convert to (ra, dec)

    // The GPS singleton has current time and orientation
    static astro::GPS* gps = fluxSvc->GPSinstance();
    double time = gps->time();

    Vector Mc_t0(McXDir, McYDir, McZDir);

    CLHEP::HepRotation R ( gps->transformToGlast(time, astro::GPS::CELESTIAL) );

    astro::SkyDir mcdir( - (R.inverse() * Mc_t0 ) );
    m_mcRa   = mcdir.ra();
    m_mcDec  = mcdir.dec();
    m_mcL = mcdir.l();
    m_mcB = mcdir.b();

    CLHEP::HepRotation Rzen ( gps->transformToGlast(time, astro::GPS::ZENITH) );
    Vector zenith = -(Rzen.inverse() * Mc_t0);
    m_mcZen  = (float) zenith.theta()*180./M_PI;
    // zero azimuth points north!
    m_mcAzim = -(float) zenith.phi()*180./M_PI + 90.0;
    if(m_mcAzim<0) m_mcAzim += 360.;

    return;
}
//------------------------------------------------------------------------------

