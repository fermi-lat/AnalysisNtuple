/** 
 @file AnalysisNtuple_load.cxx
 @brief This is needed for forcing the linker to load all components
 of the library.

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/Dll/AnalysisNtuple_load.cxx,v 1.6 2004/01/09 20:40:22 usher Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(AnalysisNtuple) {

    DECLARE_ALGORITHM( AnalysisNtupleAlg );

    DECLARE_TOOL(      TkrValsTool      );
    DECLARE_TOOL(      CalValsTool      );
    DECLARE_TOOL(      AcdValsTool      );
    DECLARE_TOOL(      McValsTool       );
    DECLARE_TOOL(      McAnalValsTool   );
    DECLARE_TOOL(      GltValsTool      );
    DECLARE_TOOL(      TkrHitValsTool   );
    DECLARE_TOOL(      VtxValsTool      );
    DECLARE_TOOL(      EvtValsTool      );
	DECLARE_TOOL(      EGRETValsTool    );

} 



