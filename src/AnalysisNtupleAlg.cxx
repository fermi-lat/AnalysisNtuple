/** @file AnalysisNtupleAlg.cxx

@brief Uses the XxxValsTools to produce a comprehensive ntuple

@author Leon Rochester



$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/AnalysisNtupleAlg.cxx,v 1.53 2012/10/03 14:12:56 bruel Exp $

*/



// Gaudi system includes

#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/AlgFactory.h"

#include "GaudiKernel/IDataProviderSvc.h" 

#include "GaudiKernel/SmartDataPtr.h"

#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/ISvcLocator.h"



// ntuple interface

#include "ntupleWriterSvc/INTupleWriterSvc.h"



// analysis tools

#include "AnalysisNtuple/IValsTool.h"



// for access to geometry perhaps

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"



#include "Event/TopLevel/EventModel.h"

#include "Event/TopLevel/Event.h"



#include <algorithm>





/** @class NtupleVisitor

@brief fields the callback from the tools; fills ntuple with names and values

@author Leon Rochester

*/

class NtupleVisitor : virtual public IValsTool::Visitor

{

public:

    /// constructor can set ntuple service and ntupleName

    NtupleVisitor(INTupleWriterSvc* ntupleSvc=0, std::string ntupleName="") 

        : m_ntupleSvc(ntupleSvc), m_ntupleName(ntupleName) {}



// LSR 14-Jul-08 code for ntuple types



    virtual IValsTool::Visitor::eVisitorRet 

        analysisValue(std::string varName, const double& value) const;

    virtual IValsTool::Visitor::eVisitorRet 

        analysisValue(std::string varName, const float& value) const;

    virtual IValsTool::Visitor::eVisitorRet 

        analysisValue(std::string varName, const int& value) const;

    virtual IValsTool::Visitor::eVisitorRet 

        analysisValue(std::string varName, const unsigned int& value) const;

    virtual IValsTool::Visitor::eVisitorRet 

        analysisValue(std::string varName, const unsigned long long& value) const;

    virtual IValsTool::Visitor::eVisitorRet 

        analysisValue(std::string varName, const char* value) const;

    virtual ~NtupleVisitor() {}

    

private:

    /// pointer to the ntuple service

    INTupleWriterSvc* m_ntupleSvc;

    /// name of the ntuple; should be the same as is set in NtupleWriterSvc

    std::string m_ntupleName;

};



// LSR 14-Jul-08 code for ntuple types



IValsTool::Visitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName,

                                                             const double& value) const

{ 

    StatusCode sc;

    if (m_ntupleSvc) {

        sc = m_ntupleSvc->addItem(m_ntupleName,  varName, &value );

        if (sc.isFailure()) return IValsTool::Visitor::ERROR;

    }    

    return IValsTool::Visitor::CONT;

}



IValsTool::Visitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName,

                                                             const float& value) const

{ 

    StatusCode sc;

    if (m_ntupleSvc) {

        sc = m_ntupleSvc->addItem(m_ntupleName,  varName, &value );

        if (sc.isFailure()) return IValsTool::Visitor::ERROR;

    }    

    return IValsTool::Visitor::CONT;

}



IValsTool::Visitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName,

                                                             const int& value) const

{ 

    StatusCode sc;

    if (m_ntupleSvc) {

        sc = m_ntupleSvc->addItem(m_ntupleName,  varName, &value );

        if (sc.isFailure()) return IValsTool::Visitor::ERROR;

    }    

    return IValsTool::Visitor::CONT;

}



IValsTool::Visitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName,

                                                             const unsigned int& value) const

{ 

    StatusCode sc;

    if (m_ntupleSvc) {

        sc = m_ntupleSvc->addItem(m_ntupleName,  varName, &value );

        if (sc.isFailure()) return IValsTool::Visitor::ERROR;

    }    

    return IValsTool::Visitor::CONT;

}





IValsTool::Visitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName,

                                                             const unsigned long long& value) const

{ 

    StatusCode sc;

    if (m_ntupleSvc) {

        sc = m_ntupleSvc->addItem(m_ntupleName,  varName, &value );

        if (sc.isFailure()) return IValsTool::Visitor::ERROR;

    }    

    return IValsTool::Visitor::CONT;

}



IValsTool::Visitor::eVisitorRet NtupleVisitor::analysisValue(std::string varName,

                                                             const char* value) const

{ 

    StatusCode sc;

    if (m_ntupleSvc) {

        sc = m_ntupleSvc->addItem(m_ntupleName,  varName, value );

        if (sc.isFailure()) return IValsTool::Visitor::ERROR;

    }    

    return IValsTool::Visitor::CONT;

}



/** @class AnalysisNtupleAlg

@brief fills the ntuple from the XxxValsTools

@author Leon Rochester

*/



class AnalysisNtupleAlg : public Algorithm {

public:

    AnalysisNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator);

    /// set parameters and attach to various perhaps useful services.

    StatusCode initialize();

    /// process one event

    StatusCode execute();

    /// clean up

    StatusCode finalize();

    

private:

    /// local utility methods

    void fixLoadOrder();

    void removeMc();

    void printHeader(MsgStream& log);



    /// number of times called

    double m_count; 

    

    /// access the ntupleWriter service to write out to ROOT ntuples

    INTupleWriterSvc *m_ntupleSvc;

    /// parameter to store the logical name of the ROOT file to write to

    std::string m_tupleName;

    

    /// Common interface to analysis tools

    std::vector<IValsTool*> m_toolvec;

    /// tool names

    std::vector<std::string> m_toolnames;

    /// switch to turn on ntupleWriterSvc, for test purposes

    bool m_doNtuple;

    bool m_doDebug;

    bool m_countCalcs;

    bool m_realData;



    // ADW 26-May-2011: Specify production tuple flag

    bool m_proTuple;



    IDataProviderSvc* m_pEventSvc;

    

    IValsTool::Visitor* m_visitor; 

};



//static const AlgFactory<AnalysisNtupleAlg>  Factory;

//const IAlgFactory& AnalysisNtupleAlgFactory = Factory;

DECLARE_ALGORITHM_FACTORY(AnalysisNtupleAlg);



AnalysisNtupleAlg::AnalysisNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator)

:Algorithm(name, pSvcLocator)

,m_count(0)

{

    // declare properties with setProperties calls

    declareProperty("tupleName",  m_tupleName="MeritTuple"); 

    // so it looks like NTupleWriterSvc property, no harm having both!

    //declareProperty("tuple_name",  m_tupleName="");   

    // List of tools to use -- maybe a bit kludgy, since the spelling needs to be correct!

    declareProperty("toolList", m_toolnames);

    declareProperty("doNtuple", m_doNtuple=true);

    declareProperty("enableDebugCalc", m_doDebug=false);

    declareProperty("countCalcs", m_countCalcs=false);  

    declareProperty("realData", m_realData=false);  



    // ADW 26-May-2011: Specify production tuple flag

    declareProperty("proTuple",m_proTuple = false);

}



StatusCode AnalysisNtupleAlg::initialize(){

    StatusCode   sc = StatusCode::SUCCESS;

    StatusCode fail = StatusCode::FAILURE;

    

    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "initialize" << endreq;



    // calc tools - default is the full set!



    // Use the Job options service to set the Algorithm's parameters

    m_toolnames.clear();



    setProperties();



    //probably a better way to do this!

    // default set: Mc is now called after Evt

    std::string toolnames [] = {"Glt", "Acd", "Acd2", "TkrHit", "Tkr", "Tree", "Vtx",  "Cal",  "Evt", "Obf", "Mc", "McTkrHit", ""};

    unsigned int i;

    unsigned int namesSize = m_toolnames.size();



    // use the default

    if(m_toolnames.empty()) {

        for (i=0; ; ++i) {

            if (toolnames[i]=="") break;

            m_toolnames.push_back(toolnames[i]);

        }

    // fresh list, use it

    } else if (m_toolnames.size()>0&&m_toolnames[0]!="+") {

        for (i=0; i<namesSize; ++i) {

            m_toolnames[i] = m_toolnames[i];

        }

    // add the input to the regular list

    } else if (namesSize>1&&m_toolnames[0]=="+") {

        m_toolnames.erase(m_toolnames.begin());

        // first the input

        namesSize = m_toolnames.size();

        for (i=0; i<namesSize; ++i) {

            m_toolnames[i] = m_toolnames[i];

        }

        // then the regulars

        for (i=0; ; ++i) {

            if (toolnames[i]=="") break;

            m_toolnames.push_back(toolnames[i]);

        }

    }



    log << MSG::INFO << endreq;

    namesSize = m_toolnames.size();

    log << MSG::INFO << namesSize << " Tools requested: ";

    for (i=0; i<namesSize; ++i) {

        log << m_toolnames[i]+"ValsTool" << " " ;

    }

    log << endreq;

    

    if( m_tupleName.empty()) {

        log << MSG::INFO << "tupleName property not set!  No ntuple output"<<endreq;

    }

    // set up tools

    IToolSvc* pToolSvc = 0;

    

    sc = service("ToolSvc", pToolSvc, true);

    if (!sc.isSuccess ()){

        log << MSG::ERROR << "Can't find ToolSvc, will quit now" << endreq;

        return StatusCode::FAILURE;

    }

    // for real data, remove all tools that start with "Mc"

    if(m_realData) removeMc();



    // take care of the Acd/Tkr load order: Acd comes first!

    // also, Mc after Evt

    // I wish we didn't have this restriction!

    fixLoadOrder();



    namesSize = m_toolnames.size();

    for (i =0; i!=namesSize; ++i){

        m_toolvec.push_back(0);

        m_toolnames[i]+="ValsTool";

        sc = pToolSvc->retrieveTool(m_toolnames[i], m_toolvec.back());

        m_toolvec.back()->setLoadOrder(i);

        if( sc.isFailure() ) {

            log << MSG::ERROR << "Unable to find tool: " 

                << m_toolnames[i] << endreq;

            return sc;

        }

    }

        

    // get a pointer to our ntupleWriterSvc

    m_ntupleSvc = 0;

    if (!m_tupleName.empty()&& m_doNtuple) {

        if (service("RootTupleSvc", m_ntupleSvc, true).isFailure()) {

            log << MSG::ERROR 

                << "AnalysisNtupleAlg failed to get the RootTupleSvc" 

                << endreq;

            return fail;

        }

    }

    

    m_visitor = new NtupleVisitor(m_ntupleSvc, m_tupleName);

    

    if (!m_tupleName.empty() && m_doNtuple) {

                

        int size = m_toolvec.size();

        for( int i =0; i<size; ++i){

            if(m_toolvec[i]->traverse(m_visitor, false, m_proTuple)==IValsTool::Visitor::ERROR) {

                log << MSG::ERROR << m_toolvec[i] << " traversal failed" << endreq;

                return fail;

            }

        }

        IDataProviderSvc* eventsvc = 0;

          sc = serviceLocator()->service( "EventDataSvc", eventsvc, true );

          if(sc.isFailure()){

              log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;

              return sc;

          }

          m_pEventSvc = eventsvc;



    }

    

    return sc;

}



StatusCode AnalysisNtupleAlg::execute()

{

    StatusCode   sc = StatusCode::SUCCESS;

    //StatusCode fail = StatusCode::FAILURE;

    

    MsgStream   log( msgSvc(), name() );



    bool countCalc = m_countCalcs;





    ++m_count;

    m_ntupleSvc->storeRowFlag(m_tupleName, true);  // needed to save the event with RootTupleSvc



    int toolCounter = 0;

    bool isException = false;

    std::vector<IValsTool*>::iterator i = m_toolvec.begin();

    for( ; i != m_toolvec.end(); ++i, ++toolCounter){

        try {

            (*i)->doCalcIfNotDone();

        } catch( std::exception& e) {

            printHeader(log);

            log << "Non-propagator exception from " << m_toolnames[toolCounter]<< endreq; 

            log << e.what() << endreq;

        } catch (...) {

            printHeader(log);

            log << "Non-propagator exception from " << m_toolnames[toolCounter]<< endreq;

            isException = true;

        }

    }

    if(isException) {



      SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);

      if (header) header->setAnalysisNtupleError();



        for(i=m_toolvec.begin() ; i != m_toolvec.end(); ++i){

            (*i)->zeroVals();

        }

        return sc;

    }





    // all the tools have been called at this point, so from now on,

    ///  we can call them with the no-calculate flag

    

    bool debugStuff = m_doDebug;



    if (countCalc || debugStuff) {

        if (m_count < 5) {

            log << MSG::INFO << "number of calculations per event: " << endreq; 

            unsigned int i;

            for(i =0; i<m_toolvec.size(); ++i){

                log << MSG::INFO << m_toolnames[i] << ": " << m_toolvec[i]->getCalcCount()<< endreq;

            }

            log << MSG::INFO << endreq;

        }else {

            log << MSG::INFO << "call count message suppressed. " << endreq;

        } 

    }



    if(debugStuff) {



        int namesSize = m_toolnames.size();

        int i;



        //do a browse



        for (i=0;i<namesSize; ++i) {

            log << MSG::INFO << "Dump of variables in " << m_toolnames[i] << endreq;

            m_toolvec[i]->browse(log);

        }



        std::string varname;

        std::vector<std::string> varnames;

        varnames.clear();

        for (i=0; i<namesSize; ++i) {

            std::string toolname = m_toolnames[i];

            if      (toolname=="McValsTool"     ) {varname = "McXErr";}

            else if (toolname=="GltValsTool"    ) {varname = "GltTotal";}

            else if (toolname=="TkrHitValsTool" ) {varname = "TkrHitsInLyr00";}

            else if (toolname=="TkrValsTool"    ) {varname = "TkrSumKalEne";}

            else if (toolname=="VtxValsTool"    ) {varname = "VtxZDir";}

            else if (toolname=="CalValsTool"    ) {varname = "CalEnergyRaw";}

            else if (toolname=="CalMipValsTool" ) {varname = "CalMipNum";}

            //else if (toolname=="GcrSelectValsTool" ) {varname = "GcrSelect[1536]","InferedZ";}

            //else if (toolname=="GcrReconValsTool" )  {varname = "GcrRecon[1536]";}

            else if (toolname=="AcdValsTool"    ) {varname = "AcdTileCount";}

            else if (toolname=="Acd2ValsTool"   ) {varname = "Acd2TileCount";}

            else if (toolname=="EvtValsTool"    ) {varname = "EvtEnergyRaw";}

            else if (toolname=="McAnalValsTool" ) {varname = "McaPrmEnegy";}

            else                                  {varname = "";}

            varnames.push_back(varname);

        }

               

        // check browse() against getVal() for each tool

                

        int vecSize = m_toolvec.size();

        for(i=0; i<vecSize; ++i){

            varname = varnames[i];

            if (varname=="") continue;

            std::string answerString;

            sc = m_toolvec[i]->getVal(varname, answerString, NOCALC);

            log  << MSG::INFO << varname << " = " << answerString << " " << endreq;  

            m_toolvec[i]->browse(log, varnames[i]);

        }        

    }     

    log << MSG::DEBUG;

    if (log.isActive()) {

        log << MSG::DEBUG << "number of calculations per event: " << endreq; 

        unsigned int i;

        for(i =0; i<m_toolvec.size(); ++i){

            log << MSG::DEBUG << m_toolnames[i] << ": " << m_toolvec[i]->getCalcCount()<< endreq;

        }

    }

    log << endreq;

    return sc;

}



StatusCode AnalysisNtupleAlg::finalize(){

    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;

    

    return sc;

}



void AnalysisNtupleAlg::fixLoadOrder()

{

    bool fixedAcd1 = false;

    bool fixedAcd2 = false;

    bool fixedMc   = false;



    MsgStream log(msgSvc(), name());

    std::vector<std::string>::iterator tkrIter, acdIter, endIter, mcIter, evtIter, tkrHitIter;



    endIter = m_toolnames.end();

    tkrIter = find(m_toolnames.begin(), m_toolnames.end(), "Tkr");

    tkrHitIter = find(m_toolnames.begin(), m_toolnames.end(), "TkrHit");

    if(tkrHitIter!=endIter&&tkrHitIter<tkrIter) tkrIter = tkrHitIter;

    acdIter = find(m_toolnames.begin(), m_toolnames.end(), "Acd");



    if(acdIter!=endIter&&tkrIter!=endIter&&acdIter-tkrIter>0) {

        m_toolnames.erase(acdIter);

        m_toolnames.insert(tkrIter, "Acd");

        fixedAcd1 = true;

    }



    endIter = m_toolnames.end();

    tkrIter = find(m_toolnames.begin(), m_toolnames.end(), "Tkr");

    tkrHitIter = find(m_toolnames.begin(), m_toolnames.end(), "TkrHit");

    if(tkrHitIter!=endIter&&tkrHitIter<tkrIter) tkrIter = tkrHitIter;

    acdIter = find(m_toolnames.begin(), m_toolnames.end(), "Acd2");



    if(acdIter!=endIter&&tkrIter!=endIter&&acdIter-tkrIter>0) {

        m_toolnames.erase(acdIter);

        m_toolnames.insert(tkrIter, "Acd2");

        fixedAcd2 = true;

    }



    mcIter  = find(m_toolnames.begin(), m_toolnames.end(), "Mc");

    evtIter = find(m_toolnames.begin(), m_toolnames.end(), "Evt");

    endIter = m_toolnames.end();



    if(mcIter!=endIter&&evtIter!=endIter&&evtIter-mcIter>0) {

        m_toolnames.push_back("Mc");

        mcIter  = find(m_toolnames.begin(), m_toolnames.end(), "Mc");

        m_toolnames.erase(mcIter);

        fixedMc = true;

    }



    if(fixedAcd1||fixedAcd2||fixedMc) {

        log << MSG::WARNING << endreq << 

            "Load order of ValsTools changed" << endreq;

        if(fixedAcd1) log << "AcdValsTool inserted before Tkr[Hit]ValsTool" << endreq;

        if(fixedAcd2) log << "Acd2ValsTool inserted before Tkr[Hit]ValsTool" << endreq;

        if(fixedMc)   log << "McValsTool moved to end" << endreq;

        unsigned int i;

        unsigned int namesSize = m_toolnames.size();

        log << MSG::WARNING << "New order: " << endreq;

        for (i=0; i<namesSize; ++i) {

            log << m_toolnames[i]+"ValsTool" << " " ;

        }

        log << endreq;

    }



}



void AnalysisNtupleAlg::removeMc() 

{

    MsgStream log(msgSvc(), name());

 //   std::vector<std::string>::iterator  endIter, beginIter, lastIter, listIter;

 //   //mckIter = find(m_toolnames.begin(), m_toolnames.end(), "McKludge");

 //   endIter   = m_toolnames.end();

        //beginIter = m_toolnames.begin();

 //   listIter = endIter;

 //   --listIter;



    unsigned int origSize = m_toolnames.size();



    unsigned int i;

        

    // one more try

    std::vector<std::string> temp;

    for(i=0;i<origSize;++i) {

        if(m_toolnames[i].find("Mc")==std::string::npos) temp.push_back(m_toolnames[i]);

    }



    m_toolnames = temp;



    unsigned int namesSize = m_toolnames.size();

    if(origSize>namesSize) {

        log << MSG::WARNING << endreq <<

            "Real Data Run: Mc tools removed" << endreq

            << "Final list (" << namesSize << " tools remain):"<< endreq;



        for (i=0; i<namesSize; ++i) {

            log << m_toolnames[i]+"ValsTool" << " " ;

        }

        log << endreq << endreq;

    }

    return;

}



void AnalysisNtupleAlg::printHeader(MsgStream& log)

{

      SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);

      unsigned long evtId = (header) ? header->event() : 0;

      long runId = (header) ? header->run() : -1;

      log << MSG::WARNING << "Caught exception (run,event): ( " 

          << runId << ", " << evtId << " ) " << endreq;

}

