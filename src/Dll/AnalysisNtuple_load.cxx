/** 
 @file AnalysisNtuple_load.cxx
 @brief This is needed for forcing the linker to load all components
 of the library.

  $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/Dll/AnalysisNtuple_load.cxx,v 1.14 2007/03/17 00:38:37 usher Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(AnalysisNtuple) {

    DECLARE_ALGORITHM( AnalysisNtupleAlg );
    DECLARE_ALGORITHM( FT1Alg            );
    DECLARE_ALGORITHM( McCoordsAlg       );
    DECLARE_ALGORITHM( ObfCoordsAlg      );

    DECLARE_TOOL(      TkrValsTool      );
    DECLARE_TOOL(      CalValsTool      );
    DECLARE_TOOL(      AcdValsTool      );
    DECLARE_TOOL(      McValsTool       );
    DECLARE_TOOL(      McAnalValsTool   );
    DECLARE_TOOL(      GltValsTool      );
    DECLARE_TOOL(      TkrHitValsTool   );
    DECLARE_TOOL(      VtxValsTool      );
    DECLARE_TOOL(      EvtValsTool      );
    //@@@FP 07/08/05
    DECLARE_TOOL(      CalMipValsTool   );
    DECLARE_TOOL(      ObfValsTool      );   //@@@LSR 03/14/07
    DECLARE_TOOL(      McTkrHitValsTool );
    //@@@FP 07/08/05
    // removed 12/19/06 LSR
    //@@@ CL 06/26/06
    //DECLARE_TOOL(      GcrSelectValsTool  );
    //DECLARE_TOOL(      GcrReconValsTool   );


} 
