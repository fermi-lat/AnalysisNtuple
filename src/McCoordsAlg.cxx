/** @file McCoordsAlg.cxx
@brief Declaration and implementation of Gaudi algorithm McCoordsAlg

$Header: /nfs/slac/g/glast/ground/cvs/merit/src/meritAlg/McCoordsAlg.cxx,v 1.16 2006/11/27 22:58:39 burnett Exp $
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
class McCoordsworker;

namespace { // anonymous namespace for file-global
    //INTupleWriterSvc* rootTupleSvc;
    IFluxSvc* fluxSvc;
    unsigned int nbOfEvtsInFile(100000);
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
    McCoordsworker * m_worker;
    //counter
    int m_count;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static const AlgFactory<McCoordsAlg>  Factory;
const IAlgFactory& McCoordsAlgFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class McCoordsworker{ 
public:
    McCoordsworker();

    void evaluate();

private:
    //std::map<std::string, double> getCelestialCoords(const Event::Exposure& exp,
    //    const CLHEP::Hep3Vector glastDir);

    bool useVertex(){ //TODO: implement
        return false;
    }
/*
    template <typename T>
        void addItem(std::string name, const T & value)
    {
        rootTupleSvc->addItem(treename, name, &value);
    }
*/
    /*
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class Item {
    public:
        Item(std::string name, char typecode=' ')
        {
            std::string type = rootTupleSvc->getItem(treename, name, m_pvalue);
            if( typecode==' ') {
                m_isFloat = type==rootType('F');
                if( !m_isFloat && type!=rootType('D')){
                    throw std::invalid_argument("McCoordsAlg: type of "+name+ " is not "+ rootType('F')+" or "+rootType('D'));
                }
            }else if( type!= rootType(typecode) ){
                throw std::invalid_argument("McCoordsAlg: type of "+name+ " is not "+ rootType(typecode));
            }
        }
        // Item behaves like a double
        operator double()const
        {
            return m_isFloat? *(float*)m_pvalue : *(double*)m_pvalue;
        }

        static std::string rootType(char code){
            if( code=='i') return "UInt_t";
            if( code=='I') return "Int_t";
            if( code=='F') return "Float_t";
            if( code=='D') return "Double_t";
            // todo: add more?
            return "unknown";
        }
        void* m_pvalue;
        bool m_isFloat;
    };
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template<typename T, char  typecode>
    class TypedItem : public Item {
    public:
        TypedItem(std::string name): Item(name, typecode){}
        T value() const{ return *static_cast<T*>(m_pvalue); }
        operator T()const{return value();}
    };

    */
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    declareProperty("NbOfEvtsInFile", nbOfEvtsInFile=100000);

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
    m_worker = new McCoordsworker();

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

McCoordsworker::McCoordsworker()
// initialize pointers to current items
: McXDir("McXDir")
, McYDir("McYDir")
, McZDir("McZDir")
{
    //now create new items 
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

    addItem( "McRa",            m_mcRa);
    addItem( "McDec",           m_mcDec);
    addItem( "McL",             m_mcL);
    addItem( "McB",             m_mcB);
    addItem( "McZenithTheta",   m_mcZen);
    addItem( "McEarthAzimuth",  m_mcAzim);
}


void McCoordsworker::evaluate()
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
    m_mcAzim = -(float) zenith.phi()*180./M_PI + 90.0;
    if(m_mcAzim<0) m_mcAzim += 360.;

    return;
}
//------------------------------------------------------------------------------

