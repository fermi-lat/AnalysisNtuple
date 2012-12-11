#include "ValBase.h"
#include <string>

static const InterfaceID IID_IPsfTool("IPsfTool", 0 , 0); 


class IPsfTool : virtual public IAlgTool

{
 public:

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IPsfTool; }

    /// load map
    virtual StatusCode loadPsf(std::string psfName) = 0;
    /// returns an interpolated value for zDir and LogE
    virtual  double computePsf(const double cl_level, 
			     const double energy,
			     const double theta, 
			     const bool isFront)  = 0;


};
