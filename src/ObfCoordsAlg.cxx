/** @file ObfCoordsAlg.cxx
@brief Declaration and implementation of Gaudi algorithm ObfCoordsAlg

$Header: /nfs/slac/g/glast/ground/cvs/merit/src/meritAlg/ObfCoordsAlg.cxx,v 1.16 2006/11/27 22:58:39 burnett Exp $
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
class ObfCoordsworker;

namespace { // anonymous namespace for file-global
    //INTupleWriterSvc* rootTupleSvc;
    IFluxSvc* fluxSvc;
    unsigned int nbOfEvtsInFile(100000);
    std::string treename("MeritTuple");
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
    ObfCoordsworker * m_worker;
    //counter
    int m_count;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static const AlgFactory<ObfCoordsAlg>  Factory;
const IAlgFactory& ObfCoordsAlgFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ObfCoordsworker{ 
public:
    ObfCoordsworker();

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
                    throw std::invalid_argument("ObfCoordsAlg: type of "+name+ " is not "+ rootType('F')+" or "+rootType('D'));
                }
            }else if( type!= rootType(typecode) ){
                throw std::invalid_argument("ObfCoordsAlg: type of "+name+ " is not "+ rootType(typecode));
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
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // tuple items expect to find
    //TypedItem<unsigned int, 'i'> EvtRun, EvtEventId;

    // these all float or double
    Item FilterXDir;
    Item FilterYDir;
    Item FilterZDir;

    //ObfCoords entries to create
    float m_obfRa,m_obfDec,m_obfL,m_obfB;
    //float m_obfZen,m_obfAzim;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ObfCoordsAlg::ObfCoordsAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)
{
    declareProperty("TreeName",  treename="MeritTuple");
    declareProperty("NbOfEvtsInFile", nbOfEvtsInFile=100000);

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
    m_worker = new ObfCoordsworker();

    service("FluxSvc", fluxSvc, true);

    m_count = 0;

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

ObfCoordsworker::ObfCoordsworker()
// initialize pointers to current items
: FilterXDir("FilterXDir")
, FilterYDir("FilterYDir")
, FilterZDir("FilterZDir")
//, FilterXhits("FilterXhits")
//, FilterYhits("FilterYhits")
{
    //now create new items 
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

    addItem( "FilterRa",            m_obfRa);
    addItem( "FilterDec",           m_obfDec);
    addItem( "FilterL",             m_obfL);
    addItem( "FilterB",             m_obfB);
    //addItem( "ObfcZenithTheta",   m_obfZen);
    //addItem( "ObfcEarthAzimuth",  m_obfAzim);
}


void ObfCoordsworker::evaluate()
{

    m_obfRa = m_obfDec = m_obfL = m_obfB = 0;

    // convert to (ra, dec)

    // The GPS singleton has current time and orientation
    static astro::GPS* gps = fluxSvc->GPSinstance();
    double time = gps->time();

    Vector filtDir(FilterXDir, FilterYDir, FilterZDir);
    if(filtDir.mag()==0) return;

    CLHEP::HepRotation R ( gps->transformToGlast(time, astro::GPS::CELESTIAL) );

    astro::SkyDir skydir( (R.inverse() * filtDir ) );
    m_obfRa   = skydir.ra();
    m_obfDec  = skydir.dec();
    m_obfL = skydir.l();
    m_obfB = skydir.b();

    //CLHEP::HepRotation Rzen ( gps->transformToGlast(time, astro::GPS::ZENITH) );
    //Vector zenith = -(Rzen.inverse() * Mc_t0);
    //m_obfZen  = (float) zenith.theta()*180./M_PI;
    //m_obfAzim = (float) zenith.phi()*180./M_PI;
    //if(m_obfAzim<0) m_obfAzim += 360.;

    return;
}
//------------------------------------------------------------------------------

