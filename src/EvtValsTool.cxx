
/** @file EvtValsTool.cxx
@brief Calculates the "Event" analysis variables from the other ntuple variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/EvtValsTool.cxx,v 1.1 2003/05/07 15:01:48 lsrea Exp $
*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "TkrUtil/ITkrGeometrySvc.h"

//#include "Event/Recon/TkrRecon/TkrVertex.h"
//#include "Event/Recon/TkrRecon/TkrFitTrack.h"
//#include "Event/Recon/AcdRecon/AcdRecon.h"

#include <algorithm>
/** @class EvtValsTool.cxx 
@brief Calculates Event Values from the other ntuple variables
*/

class EvtValsTool : public ValBase
{
public:

    EvtValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~EvtValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

private:

    //Global ACDTuple Items
    double EvtEnergyOpt;
    double EvtEnergySumOpt;
    double EvtEnergyRaw;
    double EvtMcEnergySigma;
    double EvtCalEdgeAngle;
    double EvtTkrEdgeAngle;
    double EvtLogESum;
    double EvtTkr1EFrac;
    double EvtVtxKin;
    double EvtVtxEAngle;
    double EvtTkrComptonRatio;
    double EvtTkr1EChisq;
    double EvtTkr1EFirstChisq;
    double EvtTkr1EQual;
    double EvtTkr2EChisq;
    double EvtTkr2EFirstChisq;
    double EvtTkr2EQual;

    IValsTool* m_pMcTool;
    IValsTool* m_pGltTool;
    IValsTool* m_pTkrHitTool;
    IValsTool* m_pTkrTool;
    IValsTool* m_pVtxTool;
    IValsTool* m_pCalTool;
    IValsTool* m_pAcdTool;

    ITkrGeometrySvc* pTkrGeoSvc;

};

// Static factory for instantiation of algtool objects
static ToolFactory<EvtValsTool> s_factory;
const IToolFactory& EvtValsToolFactory = s_factory;

// Standard Constructor
EvtValsTool::EvtValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this);
}

StatusCode EvtValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the services

    // try using one tool from another

    if(service( "TkrGeometrySvc", pTkrGeoSvc, true ).isFailure()) {
        log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
        return fail;
    }

    IToolSvc* pToolSvc = 0; 
    sc = service("ToolSvc", pToolSvc, true);
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Can't find ToolSvc, will quit now" << endreq;
        return StatusCode::FAILURE;
    }

    m_pMcTool = 0;
    sc = pToolSvc->retrieveTool("McValsTool", m_pMcTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "McValsTool" << endreq;
        return sc;
    }
    m_pGltTool = 0;
    sc = pToolSvc->retrieveTool("GltValsTool", m_pGltTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "GltValsTool" << endreq;
        return sc;
    }
    m_pTkrHitTool = 0;
    sc = pToolSvc->retrieveTool("TkrHitValsTool", m_pTkrHitTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "TkrHitValsTool" << endreq;
        return sc;
    }
    m_pTkrTool = 0;
    sc = pToolSvc->retrieveTool("TkrValsTool", m_pTkrTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "TkrValsTool" << endreq;
        return sc;
    }
    m_pVtxTool = 0;
    sc = pToolSvc->retrieveTool("VtxValsTool", m_pVtxTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "VtxValsTool" << endreq;
        return sc;
    }
    m_pCalTool = 0;
    sc = pToolSvc->retrieveTool("CalValsTool", m_pCalTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "CalValsTool" << endreq;
        return sc;
    }
    m_pAcdTool = 0;
    sc = pToolSvc->retrieveTool("AcdValsTool", m_pAcdTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "AcdValsTool" << endreq;
        return sc;
    }


    // load up the map

    addItem("EvtEnergyOpt",     &EvtEnergyOpt);
    addItem("EvtEnergySumOpt",  &EvtEnergySumOpt);
    addItem("EvtEnergyRaw",     &EvtEnergyRaw);
    addItem("EvtMcEnergySigma", &EvtMcEnergySigma);
    addItem("EvtCalEdgeAngle",  &EvtCalEdgeAngle);
    addItem("EvtTkrEdgeAngle",  &EvtTkrEdgeAngle);
    addItem("EvtLogESum",       &EvtLogESum);
    addItem("EvtTkr1EFrac",     &EvtTkr1EFrac);
    addItem("EvtVtxKin",        &EvtVtxKin);
    addItem("EvtVtxEAngle",     &EvtVtxEAngle);
    addItem("EvtTkrComptonRatio", &EvtTkrComptonRatio);
    addItem("EvtTkr1EChisq",    &EvtTkr1EChisq);
    addItem("EvtTkr1EFirstChisq", &EvtTkr1EFirstChisq);
    addItem("EvtTkr1EQual",     &EvtTkr1EQual);
    addItem("EvtTkr2EChisq",    &EvtTkr2EChisq);
    addItem("EvtTkr2EFirstChisq", &EvtTkr2EFirstChisq);
    addItem("EvtTkr2EQual",     &EvtTkr2EQual);

    zeroVals();

    return sc;
}

StatusCode EvtValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    //SmartDataPtr<Event::AcdRecon>           pACD(m_pEventSvc,EventModel::AcdRecon::Event);

    //Make sure we have valid ACD data
    /*
    if (pACD)
    {
    ACD_Total_Energy  = pACD->getEnergy();
    } else {
    return StatusCode::FAILURE;
    }
    */


    double eTkr, eCal;
    if (m_pTkrTool->getVal("TkrEnergyCorr", eTkr).isSuccess() && m_pCalTool->getVal("CalEnergyCorr", eCal).isSuccess()) {
        EvtEnergyOpt = eTkr + eCal;
    }

    double eCalSumCorr;
    if(m_pCalTool->getVal("CalEneSumCorr", eCalSumCorr).isSuccess()) {
        EvtEnergySumOpt = eTkr + eCalSumCorr;
    }

    double eCalSum;
    if(m_pCalTool->getVal("CalEnergySum", eCalSum).isSuccess()) {
        EvtEnergyRaw = eTkr + eCalSum;
    }

    
    double mcEnergy;
    if(m_pMcTool->getVal("McEnergy", mcEnergy).isSuccess()){
        EvtMcEnergySigma = (EvtEnergySumOpt - mcEnergy)/(.1*mcEnergy);
    }

    double tkrEdge, calEdge, tkr1ZDir;
    if(m_pTkrTool->getVal("Tkr1ZDir",tkr1ZDir).isSuccess()) {
        double sTkr = sqrt(1.-tkr1ZDir*tkr1ZDir);
        if (m_pTkrTool->getVal("TkrTwrEdge", tkrEdge).isSuccess()) {
            EvtTkrEdgeAngle = (30.-tkrEdge)/sTkr;
        }
        if (m_pCalTool->getVal("CalTwrEdge", calEdge).isSuccess()) {
            EvtCalEdgeAngle = (30. -calEdge)/sTkr;
        }
    }
    
    EvtLogESum = log10(std::min(std::max(EvtEnergySumOpt,40.),15000.));

    double tkr1ConE;
    if (m_pTkrTool->getVal("Tkr1ConEne",tkr1ConE).isSuccess()) {
        EvtTkr1EFrac = tkr1ConE/EvtEnergySumOpt;
    }

    double vtxAngle;
    if (m_pVtxTool->getVal("VtxAngle", vtxAngle).isSuccess()) {
        EvtVtxKin = vtxAngle*EvtEnergySumOpt*EvtEnergySumOpt/tkr1ConE;
    }

    EvtVtxEAngle = vtxAngle*EvtEnergySumOpt;
    
    double totHits, tkr1First;
    if (m_pTkrTool->getVal("TkrTotalHits", totHits).isSuccess()) {
        if (m_pTkrTool->getVal("Tkr1FirstLayer", tkr1First).isSuccess()){
            EvtTkrComptonRatio = totHits/(2.*(pTkrGeoSvc->numLayers()-tkr1First));
        }
    }
    double logE = EvtLogESum;
    double logE2 = logE*logE;
    double tkr1Chisq;
    if (m_pTkrTool->getVal("Tkr1Chisq", tkr1Chisq).isSuccess()) {
        EvtTkr1EChisq = tkr1Chisq/(9.22 - 3.7*logE + .412*logE2);
    }
    double tkr1_1stChisq;
    if (m_pTkrTool->getVal("Tkr1FirstChisq", tkr1_1stChisq).isSuccess()) {
        EvtTkr1EFirstChisq = tkr1_1stChisq/(4.51 - 1.8*logE + .23*logE2);
    }
    double tkr1Qual;
    if (m_pTkrTool->getVal("Tkr1Qual", tkr1Qual).isSuccess()) {
        EvtTkr1EQual = tkr1Qual/(58. + 8.39*sqrt(logE-1.6)
            - 2.46*logE);
    }
    double tkr2Chisq;
    if (m_pTkrTool->getVal("Tkr2Chisq", tkr2Chisq).isSuccess()) {
        EvtTkr2EChisq = tkr2Chisq/(8.95 - 4.77*logE + .948*logE2);
    }
    double tkr2_1stChisq;
    if (m_pTkrTool->getVal("Tkr2FirstChisq", tkr2_1stChisq).isSuccess()) {
        EvtTkr2EFirstChisq = tkr2_1stChisq/(6.3 - 4.09*logE + .971*logE2);
    }
    double tkr2Qual;
    if (m_pTkrTool->getVal("Tkr2Qual", tkr2Qual).isSuccess()) {
        EvtTkr2EQual = tkr2Qual/(61.3 + 34.8*sqrt(logE-1.6)
            - 16.4*logE);
    }

    return sc;
}
