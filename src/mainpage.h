// Mainpage for doxygen

/** @mainpage package AnalysisNtuple

 @author Leon Rochester
 
 @section description Description

  The package contains a set of tools and an algorithm that runs them
  to produce a comprehensive ntuple of recon results. The tools,
  which can be invoked independently from the algorithm, also
  allow access to all of the calculated variables.

  The package name is a bit misleading, since the only the algorithm depends
  on NTupleWriterSvc, not the tools.

 @section AnalysisNtupleAlg AnalysisNtupleAlg

  This algorithm writes out a comprehensive ntuple of derived analysis quantities.
  it instantiates a set of tools, each dealing with a catagory of variables, 
  and sets up and calls the traverse() method of each tool with a pointer to a 
  visitor. Then it fills the ntuple with information from the visitor callback.

  In debug mode, it exercises some other methods of the tools: getVal(), to return
  an individual value, and browse(), to print out a name and value.

  Since all the tools have the same abstract interface, the manipulation of these
  tools can be handled in a loop.

  AnalysisNtupleAlg also serves as the test routine for this package.

 @section tools Tools

  There are currently seven tools: McValsTool, GltValsTool, CalValsTool, 
  TkrValsTool, TkrHitValsTool, VtxValsTool, CalValsTool, AcdValsTool and EvtValsTool.

  EvtValsTool is special, in that it produces ntuple variables from the variables in 
  the other tools, using the getVal() method of those tools. These are generally needed 
  for the Gleam implementation of the classification trees.

  The common interface is IValsTool, which provides access to the values in each tool,
  and handles the interactions with the visitor through the class ValsVisitor.

  Each tool fills a list of variable names and pointers to the values they represent. 
  This is what allows the visitor to return names and values in a loop.

  Another method of each tool receives a call from the Gaudi IncidentSvc, to signal
  the beginning of a new event. This signal allows the tool to do its calculations 
  at most once for each event.

 @section jobOptions jobOptions
  
 @param AnalysisNtupleAlg.tupleName 
  Sets the name of the ntuple. (Default is "".) Must be the one opened by NTupleWriteSvc,
  but there's no way for the code to find out what that is.  The name for
  a standard GLEAM job is "GLEAM".
 
 @param AnalysisNtupleAlg.tuple_name
  An alternate spelling, to match the one used by NTupleWriterSvc.

 @param AnalysisNtupleAlg.toolList
  A vector of strings that sets the list of tools requested. (Default is the 8 standard tools.) 
  Each string is the name
  of a tool with the "ValsTool" removed, so "TkrHitValsTool" becomes "TkrHit". This means that
  any tool written for the package should be named "BlahBlahValsTool". 
  
  <hr>
 @section notes release.notes
  release.notes
 <hr>
 @section requirements requirements
 @verbinclude requirements
 <hr>

 */

