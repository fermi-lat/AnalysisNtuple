/** 
* @file GlastDigi_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/Dll/AnalysisNtuple_load.cxx,v 1.1.1.1 2003/02/23 21:31:10 burnett Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(AnalysisNtuple) {
    DECLARE_TOOL(    TkrValsTool      );
    DECLARE_TOOL(    CalValsTool      );
    DECLARE_TOOL(    AcdValsTool      );
    DECLARE_TOOL(    MCValsTool       );
    DECLARE_TOOL(    GltValsTool      );
    DECLARE_TOOL(    TkrHitValsTool   );


} 



