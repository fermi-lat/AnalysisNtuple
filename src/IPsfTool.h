#include "ValBase.h"
#include <string>

static const InterfaceID IID_IPsfTool("IPsfTool", 0 , 0); 


class IPsfTool : virtual public IAlgTool

{
 public:

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IPsfTool; }

    /// load map
    virtual StatusCode loadPsf(const std::string psfVersion, 
			 const std::string psfPath) = 0;
    /// returns an interpolated value for zDir and LogE
    virtual  double computePsf(const double cl_level, 
			     const double energy,
			     const double theta, 
			     const bool isFront)  = 0;

    void setEstimate(const double value){m_psf_estimate=value;}

 protected:
    ///initial estimate for the Psf to be computed
    double m_psf_estimate;
};
