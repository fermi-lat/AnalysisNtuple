/** @file ObfCoordsAlg.cxx
@brief Declaration and implementation of Gaudi algorithm ObfCoordsAlg

$Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/AnalysisNtuple/src/ObfCoordsAlg.cxx,v 1.16.80.1 2011/04/19 18:10:25 heather Exp $
*/
// Include files

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "astro/GPS.h"
#include "geometry/Vector.h"

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

//static const AlgFactory<ObfCoordsAlg>  Factory;
//const IAlgFactory& ObfCoordsAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(ObfCoordsAlg);

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

    // get the GPS instance
    gps = astro::GPS::instance();

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
<tr><td> GrbRa, GrbDec 
<td>F<td>  (deg) reconstructed direction in equatorial coordinates       
<tr><td> GrbL, GrbB 
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
    addItem( "GrbRa",     m_grbRa);
    addItem( "GrbDec",    m_grbDec);
    addItem( "GrbL",      m_grbL);
    addItem( "GrbB",      m_grbB);
}


void ObfCworker::evaluate()
{
    m_grbRa = m_grbDec = m_grbL = m_grbB = 0;

    Vector grbDir(GrbXDir, GrbYDir, GrbZDir);
    if (grbDir.mag()==0) return;
    astro::SkyDir skyGrbdir( gps->toSky(-grbDir) );
    m_grbRa  = skyGrbdir.ra();
    m_grbDec = skyGrbdir.dec();
    m_grbL   = skyGrbdir.l();
    m_grbB   = skyGrbdir.b();

    return;
}
//------------------------------------------------------------------------------

