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

  This algorithm instantiates a set of tools, each dealing with a category of variables, 
  and sets up and calls the traverse() method of each tool with a pointer to a 
  visitor. Then it fills the ntuple with information from the visitor callback. 
  The code can handle ints, uints, floats and doubles.

  Currently, AnalysisNtupleAlg outputs the ntuple variables to the merit ntuple. It also serves as the test   routine for the package

  The definition of the variables can be found in the "Related Pages" of the doxygen output.

  In debug mode, it exercises some other methods of the tools: getVal(), to return
  an individual value, and browse(), to print out a name and value.

  Since all the tools have the same abstract interface, the manipulation of these
  tools can be handled in a loop.

 @section tools Tools

  Each tool adds variables to the ntuple. To add a new variable:

 @verbatim
  float myVar; // or int, unsigned, or double
 @endverbatim
  
  At init time:

 @verbatim
    addItem("MyVarName", &myVar);
 @endverbatim

  At execution time:
  
 @verbatim
    myVar = sin(x)/x;
 @endverbatim

  Any item in the ntuple can be retrieved:

 @verbatim
    getVal("MyVarName", &retrievedVar, NOCALC);
 @endverbatim

  NOCALC says don't recalculate the ntuple contents. The default is CALC.

  It's now also possible to define arrays:

 @verbatim
  
  double arrayVar[3];
  ...
  addItem("ArrayVar[3]", arrayVar);  // note: ArrayVar is a pointer
  ...
  for(i=0;i<3;++i) { arrayVar[i] = 3.14159*i; }
  ...
  getVal("ArrayVar[1]", &arrayVar[1], NOCALC); // 
 @endverbatim
 
  Note: because Insightful Miner (C) doesn't recognize arrays, we aren't using this feature in the
  official merit ntuple. But it's available for private variables.

  
  There are currently ten tools. Eight are invoked to produce the standard merit ntuple:
  McValsTool, GltValsTool, CalValsTool, CalMipValsTool,
  TkrValsTool,  VtxValsTool, AcdValsTool and EvtValsTool.

  Two other tools are not part of the standard job, but may be invoked by adding their 
  names to the tools list: 
  TkrHitValsTool gives some addtional information about cluster counts, 
  and McAnalValsTool does an analysis of the simulated event; it requires 
  that various nonstandard packages have been run.
  
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

 @param AnalysisNtupleAlg.doNtuple
  Controls the generation of the ntuple. (Default is true.) This package now does the output
  for the ntuple.
  
 @param AnalysisNtupleAlg.tupleName 
  Sets the name of the ntuple. (Default is "MeritTuple".) Must be the one opened by NTupleWriteSvc,
  but there's no way for the code to find out what that is.  
 
 @param AnalysisNtupleAlg.toolList
  A vector of strings that sets the list of tools requested. (Default is the 8 standard tools.) 
  Each string is the name
  of a tool with the "ValsTool" removed, so "TkrHitValsTool" becomes "TkrHit". This means that
  any tool written for the package should be named "BlahBlahValsTool". 

 @param AnalysisNtupleAlg.enableDebugCalc
  Does some tests on the ntuple vars for every event... doesn't seem to be working too well at
  the moment!

 
 @param AnalysisNtupleAlg.countCalcs
  Counts the number of times that a given Tool is called for each event
  
 @param PtValsAlg.pointing_info_tree_name ["MeritTuple"] Sets the name of the tuple to write Pt info to
 @param PtValsAlg.PointingHistory [""]
  name of a FT2 FITS file, or a pointing history ascii file. If not set, and FluxSvc
  is not used, the default orbit will be used to define the pointing

  @param FT1Alg.EnableAberrationCorrection [false]
  if set true, the flag in GPS will be set to apply the orbital aberration correction when converting from 
  instrument to celestical coordinates.
  <hr>
 @section vars Description of the variables
 Here are the <A HREF="./anatup_vars.html"> standard</A> variables, and here are the 
 <A HREF="./anatup_vars_optional.html"> optional</A> ones.

  <hr>
 @section notes Release Notes
  release.notes
 <hr>
 @section requirements requirements
 @verbinclude requirements
 <hr>

 */
