/** @file McCoordsAlg.cxx
@brief Declaration and implementation of Gaudi algorithm McCoordsAlg

$Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/AnalysisNtuple/src/McCoordsAlg.cxx,v 1.11 2011/12/12 20:36:35 heather Exp $
*/
// Include files

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "astro/GPS.h"
#include "geometry/Vector.h"

#include <cassert>
#include <map>
#include <stdexcept>

// forward declatation of the worker
class McCworker;

namespace { // anonymous namespace for file-global
    std::string treename("MeritTuple");
    astro::GPS* gps;
    const double R2D = 180./M_PI;

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
    // flag for no MC data, skip this alg
    bool m_noMC;

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//static const AlgFactory<McCoordsAlg>  Factory;
//const IAlgFactory& McCoordsAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(McCoordsAlg);

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

	m_noMC = false;

    // get a pointer to RootTupleSvc 
    if( (sc = service("RootTupleSvc", rootTupleSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
        return sc;
    }
	try {
        m_worker = new McCworker();
	}
	catch(std::exception ex) {
		m_noMC = true;
                log << MSG::INFO << "There are no Mc[X/Y/Z]Dir variables, so McCoordsAlg will be skipped" << endreq;
		return sc;
	}

    // get the GPS instance
    gps = astro::GPS::instance();

    m_count = 0;

    return sc;
}

StatusCode McCoordsAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    m_count++;

	if(m_noMC) return sc;

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

    Vector Mc_t0(McXDir, McYDir, McZDir);
    // Mc particle direction points down... 
    // toSky converts a *particle* direction
    // into a direction on the sky, so the minus-sign is taken care of!

    astro::SkyDir mcdir = gps->toSky( Mc_t0 );
    m_mcRa   = mcdir.ra();
    m_mcDec  = mcdir.dec();
    m_mcL = mcdir.l();
    m_mcB = mcdir.b();

    // Local zenith coordinates
    astro::SkyDir zenith(gps->zenithDir());  // pointing direction
    double zenith_theta = mcdir.difference(zenith); 
    if( fabs(zenith_theta)<1e-8) zenith_theta=0;
    // all this to do the azimuth angle :-(
    CLHEP::Hep3Vector north_pole(0,0,1);
    CLHEP::Hep3Vector east_dir( north_pole.cross(zenith()).unit() ); // east is perp to north_pole and zenith
    CLHEP::Hep3Vector north_dir( zenith().cross(east_dir));
    double earth_azimuth=atan2( mcdir().dot(east_dir), mcdir().dot(north_dir) );
    if( earth_azimuth <0) earth_azimuth += 2*M_PI; // to 0-360 deg.
    if( fabs(earth_azimuth)<1e-8) earth_azimuth=0;
    m_mcZen  = zenith_theta*R2D;;
    m_mcAzim = earth_azimuth*R2D;

    return;
}
//------------------------------------------------------------------------------

