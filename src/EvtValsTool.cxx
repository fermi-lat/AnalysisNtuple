
/** @file EvtValsTool.cxx
@brief Calculates the "Event" analysis variables from the other ntuple variables
@author Bill Atwood, Leon Rochester

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/EvtValsTool.cxx,v 1.11 2003/11/13 20:32:25 lsrea Exp $
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
    //double EvtEnergyOpt;
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
    double EvtTkrEComptonRatio;
    double EvtPSFModel; 

    double EvtTkr1EChisq;
    double EvtTkr1EFirstChisq;
    double EvtTkr1EQual;
    double EvtTkr1ECovDet; 
    double EvtTkr1PSFMdRat;

    double EvtTkr2EChisq;
    double EvtTkr2EFirstChisq;
    double EvtTkr2EQual;

    double EvtCalETLRatio;
    double EvtCalEXtalRatio;
    double EvtCalEXtalTrunc;
    double EvtCalETrackDoca;
    double EvtCalETrackSep;
    double EvtVtxEEAngle;
    double EvtVtxEDoca;
    double EvtVtxEHeadSep;

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

    //addItem("EvtEnergyOpt",     &EvtEnergyOpt);
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
    addItem("EvtTkrEComptonRatio", &EvtTkrEComptonRatio);
    addItem("EvtPSFModel",      &EvtPSFModel);

    addItem("EvtTkr1EChisq",    &EvtTkr1EChisq);
    addItem("EvtTkr1EFirstChisq", &EvtTkr1EFirstChisq);
    addItem("EvtTkr1EQual",     &EvtTkr1EQual);
    addItem("EvtTkr1PSFMdRat",  &EvtTkr1PSFMdRat);
    addItem("EvtTkr1ECovDet",   &EvtTkr1ECovDet);

    addItem("EvtTkr2EChisq",    &EvtTkr2EChisq);
    addItem("EvtTkr2EFirstChisq", &EvtTkr2EFirstChisq);
    addItem("EvtTkr2EQual",     &EvtTkr2EQual);

    addItem("EvtCalETLRatio",   &EvtCalETLRatio);
    addItem("EvtCalEXtalRatio", &EvtCalEXtalRatio);
    addItem("EvtCalEXtalTrunc", &EvtCalEXtalTrunc);
    addItem("EvtCalETrackDoca", &EvtCalETrackDoca);
    addItem("EvtCalETrackSep",  &EvtCalETrackSep);
    addItem("EvtVtxEEAngle",    &EvtVtxEEAngle);
    addItem("EvtVtxEDoca",      &EvtVtxEDoca);
    addItem("EvtVtxEHeadSep",   &EvtVtxEHeadSep);

    zeroVals();

    return sc;
}

StatusCode EvtValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // we may want to add TDS stuff to this method, but we haven't needed it yet.

    int firstCheck = m_check;
    int nextCheck = -1;

    // since we know what's happening here, we can plan a little
    // the idea is to call the first check of each tool with the called check value,
    // and the rest with the no-calc value
    // so be careful when adding calls or moving stuff around!!!!

    double eCalSum, eTkr;
    if(    m_pCalTool->getVal("CalEnergySum", eCalSum, firstCheck).isSuccess()
        && m_pTkrTool->getVal("TkrEnergyCorr", eTkr, firstCheck).isSuccess()) {
        EvtEnergyRaw = eTkr + eCalSum;
    }

    double eTkrKalEne, eCalRLn, eTkrBest; //, eCal
    int CAL_Type;
    if (  // m_pCalTool->getVal("CalEnergyCorr", eCal, firstCheck).isSuccess() &&
        m_pCalTool->getVal("CalTotRLn", eCalRLn, firstCheck).isSuccess()
        && m_pTkrTool->getVal("TkrSumKalEne", eTkrKalEne, firstCheck).isSuccess()) 
    {
        eTkrBest = std::max(eTkr+eCalSum, eTkrKalEne);
        if(eCalSum < 100 || eCalRLn < 2) {
            if(eCalSum < 5 || eCalRLn < 2) {
                CAL_Type = 0;
            }
            else {CAL_Type = 1;}
        }
        else {CAL_Type = 2;} 

       // if(CAL_Type == 0) EvtEnergyOpt = eTkrBest;
       //   else              EvtEnergyOpt = eTkr + eCal;
    }
   
    double eCalSumCorr;
    if(m_pCalTool->getVal("CalEneSumCorr", eCalSumCorr, nextCheck).isSuccess()) {
        if(CAL_Type == 0) EvtEnergySumOpt = eTkrBest;
        else              EvtEnergySumOpt = eTkr + eCalSumCorr;
    }
    
    double mcEnergy;
    if(m_pMcTool->getVal("McEnergy", mcEnergy, firstCheck).isSuccess()){
        EvtMcEnergySigma = (EvtEnergySumOpt - mcEnergy)/(mcEnergy);
    }

    double tkrEdge, calEdge, tkr1ZDir = -1., tkr1ZDir2 = 1.; ;
    if(m_pTkrTool->getVal("Tkr1ZDir",tkr1ZDir, nextCheck).isSuccess()) {
        tkr1ZDir2 = tkr1ZDir*tkr1ZDir; 
        double sTkr = sqrt(1.-tkr1ZDir2);
        if (m_pTkrTool->getVal("TkrTwrEdge", tkrEdge, nextCheck).isSuccess()) {
            EvtTkrEdgeAngle = (30.-tkrEdge)/sTkr;
        }
        if (m_pCalTool->getVal("CalTwrEdge", calEdge, nextCheck).isSuccess()) {
            EvtCalEdgeAngle = (30. -calEdge)/sTkr;
        }
    }
    EvtPSFModel = sqrt(pow((.061/pow((std::max(EvtEnergySumOpt,1.)/100),.8)),2) + (.001745*.001745));
    EvtLogESum = log10(std::min(std::max(EvtEnergySumOpt,10.),1000000.));
    double logE = std::min(std::max(EvtLogESum,1.3), 4.7);
    double logE2 = logE*logE; 
    
    double tkr1CovDet;
    if (m_pTkrTool->getVal("Tkr1CovDet",tkr1CovDet, nextCheck).isSuccess()) {
        EvtTkr1ECovDet = tkr1CovDet/pow(std::max(EvtEnergySumOpt,1.0), 1.3);
    }

    double tkr1ThetaErr, tkr1PhiErr;
    if (m_pTkrTool->getVal("Tkr1ThetaErr",tkr1ThetaErr, nextCheck).isSuccess() &&
        m_pTkrTool->getVal("Tkr1PhiErr",tkr1PhiErr, nextCheck).isSuccess()) {
        EvtTkr1PSFMdRat = sqrt(tkr1ThetaErr*tkr1ThetaErr + tkr1PhiErr*tkr1PhiErr)/ 
                          EvtPSFModel;
    }

    double tkr1ConE;
    if (m_pTkrTool->getVal("Tkr1ConEne",tkr1ConE, nextCheck).isSuccess()) {
        if(EvtEnergySumOpt>0.0) EvtTkr1EFrac = tkr1ConE/EvtEnergySumOpt;
    }

    double vtxAngle;
    if (m_pVtxTool->getVal("VtxAngle", vtxAngle, firstCheck).isSuccess()) {
        if (tkr1ConE>0.0) EvtVtxKin = vtxAngle*EvtEnergySumOpt*EvtEnergySumOpt/tkr1ConE;
    }

    EvtVtxEAngle = vtxAngle*EvtEnergySumOpt;
    
    double totHits, tkr1First;
    if (m_pTkrTool->getVal("TkrTotalHits", totHits, nextCheck).isSuccess()) {
        if (m_pTkrTool->getVal("Tkr1FirstLayer", tkr1First, nextCheck).isSuccess()){
            EvtTkrComptonRatio = totHits/(2.*(pTkrGeoSvc->numLayers()-tkr1First));
            EvtTkrEComptonRatio = EvtTkrComptonRatio/(-3.77 + 3.57*logE - .547*logE2)
                                  /(1.69 + 1.74*tkr1ZDir + .987*tkr1ZDir2);
        }
    }

    double tkr1Chisq;
    if (m_pTkrTool->getVal("Tkr1Chisq", tkr1Chisq, nextCheck).isSuccess()) {
        EvtTkr1EChisq = tkr1Chisq/(11.5 - 3.73*logE + .362*logE2)
                                 /(.558 - 2.53*tkr1ZDir - 2.50*tkr1ZDir2);
    }
    double tkr1_1stChisq;
    if (m_pTkrTool->getVal("Tkr1FirstChisq", tkr1_1stChisq, nextCheck).isSuccess()) {
        EvtTkr1EFirstChisq = tkr1_1stChisq/(6.37 - 1.87*logE + .199*logE2)
                                /(1.99 + 1.39*tkr1ZDir + .0805*tkr1ZDir2);
    }
    double tkr1Qual;
    if (m_pTkrTool->getVal("Tkr1Qual", tkr1Qual, nextCheck).isSuccess()) {
        EvtTkr1EQual = tkr1Qual/(49.3 + 4.52*logE - .518*logE2)
                                /(1.04 + .197*tkr1ZDir + .180*tkr1ZDir2);
    }
    double tkr2Chisq;
    if (m_pTkrTool->getVal("Tkr2Chisq", tkr2Chisq, nextCheck).isSuccess()) {
        EvtTkr2EChisq = tkr2Chisq/(4.74- 1.06*logE + .257*logE2)
                                /(.366 - 3.68*tkr1ZDir - 3.18*tkr1ZDir2);
    }
    double tkr2_1stChisq;
    if (m_pTkrTool->getVal("Tkr2FirstChisq", tkr2_1stChisq, nextCheck).isSuccess()) {
        EvtTkr2EFirstChisq = tkr2_1stChisq/ (.472 + .77*logE + .0646*logE2)
                                /(.0404 - 4.57*tkr1ZDir - 3.66*tkr1ZDir2);
    }
    double tkr2Qual;
    if (m_pTkrTool->getVal("Tkr2Qual", tkr2Qual, nextCheck).isSuccess()) {
        double minLogE = std::min(logE, 2.77);
        double minLogE2 = minLogE*minLogE; 
        EvtTkr2EQual = tkr2Qual/ (-6.53 +42.1*minLogE - 7.63*minLogE2)
                                /(1.18 + .541*tkr1ZDir + .368*tkr1ZDir2);
    }

    double calTrms, calLrms;
    if(m_pCalTool->getVal("CalTransRms", calTrms, nextCheck).isSuccess()&
       m_pCalTool->getVal("CalLongRms", calLrms, nextCheck).isSuccess()){
           if(calLrms > 0.) {
             EvtCalETLRatio = (calTrms/calLrms)/ (9 - 4.01*logE + .481*logE2)
                                /(1.47 + .907*tkr1ZDir + .338*tkr1ZDir2);
           }
    }

    double calXtalRatio;
    if(m_pCalTool->getVal("CalXtalRatio", calXtalRatio, nextCheck).isSuccess()) {
        EvtCalEXtalRatio = calXtalRatio/(3.47-1.40*logE  + .146*logE2)
                               /(.960 + .0933*tkr1ZDir + .167*tkr1ZDir2);
    }

    double calXtalsTrunc;
    if(m_pCalTool->getVal("CalXtalsTrunc", calXtalsTrunc, nextCheck).isSuccess()) {
        double term1 = std::max(1., (-85.4 + 65.2*logE - 10.4*logE2));
        double logE33 = logE-3.3; 
        double term2 = (logE33 > 0.) ? std::min(14., 12.*logE33*logE33): 0.;
        EvtCalEXtalTrunc = calXtalsTrunc/(term1 + term2)/(.935 - .382*tkr1ZDir - .343*tkr1ZDir2);
    }

    double calTrackDoca;
    if(m_pCalTool->getVal("CalTrackDoca", calTrackDoca, nextCheck).isSuccess()) {
        EvtCalETrackDoca = calTrackDoca*sqrt(EvtEnergySumOpt)/(3100 - 1500*logE + 239*logE2)
                                /(5.06 + 10.5*tkr1ZDir +6.17*tkr1ZDir2);
    }

    double calTrackSep;
    if(m_pCalTool->getVal("CalTrackSep", calTrackSep, nextCheck).isSuccess()) {
        double logEmin = std::min(logE, 3.6); 
        EvtCalETrackSep = calTrackSep/(1380 - 692*logEmin + 91.3*logEmin*logEmin)
                                /(3.85 + 5.98*tkr1ZDir +2.53*tkr1ZDir2);
    }

    EvtVtxEEAngle = EvtVtxEAngle/ (96.0 -  70.7*logE + 17.6*logE2) 
                                / (1.54 + 1.46*tkr1ZDir + .889*tkr1ZDir2);
    double vtxDoca;
    if (m_pVtxTool->getVal("VtxDOCA", vtxDoca, firstCheck).isSuccess()) {
        EvtVtxEDoca = vtxDoca/(1.55 - .685*logE+ .0851*logE2) 
                                / (2.21 + 3.01*tkr1ZDir + 1.59*tkr1ZDir2);
    }

    double vtxHeadSep;
    if (m_pVtxTool->getVal("VtxHeadSep", vtxHeadSep, firstCheck).isSuccess()) {
        EvtVtxEHeadSep = vtxHeadSep/ (2.83 - .94*logE + .108*logE2) 
                                / (2.45 + 3.34*tkr1ZDir + 1.58*tkr1ZDir2);
    }

    return sc;
}
