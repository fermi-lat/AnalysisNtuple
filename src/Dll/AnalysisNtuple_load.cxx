/** 
 @file AnalysisNtuple_load.cxx
 @brief This is needed for forcing the linker to load all components
 of the library.

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/Dll/AnalysisNtuple_load.cxx,v 1.4 2003/03/05 02:03:29 lsrea Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(AnalysisNtuple) {

    DECLARE_ALGORITHM( AnalysisNtupleAlg );

    DECLARE_TOOL(      TkrValsTool      );
    DECLARE_TOOL(      CalValsTool      );
    DECLARE_TOOL(      AcdValsTool      );
    DECLARE_TOOL(      McValsTool       );
    DECLARE_TOOL(      GltValsTool      );
    DECLARE_TOOL(      TkrHitValsTool   );
    DECLARE_TOOL(      VtxValsTool      );
    DECLARE_TOOL(      EvtValsTool      );

} 



