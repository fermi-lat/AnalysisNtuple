/** @file EvtValsTool.cxx
@brief Calculates an expected PSF containment radius for each evt, based on a given psf version
@author Johann Cohen-Tanugi

$Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/PsfValsTool.cxx,v 1.1 2012/12/11 15:24:54 cohen Exp $
*/
#include <sstream>
#include <stdexcept>

#include "IPsfTool.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "facilities/Util.h"                // for expandEnvVar

#include "st_facilities/Bilinear.h"
#include "st_facilities/GaussianQuadrature.h"
#include "st_facilities/FitsTable.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

namespace {
   double sqr(double x) {
      return x*x;
   }
   class Array {
   public:
      Array(const std::vector<double> & values, size_t nx) 
         : m_values(values), m_nx(nx) {}
      double operator()(size_t iy, size_t ix) const {
         return m_values[iy*m_nx + ix];
      }
   private:
      const std::vector<double> & m_values;
      size_t m_nx;
   };
}

////////////////////////////////////////////////
class PsfObject{

public:

  PsfObject(std::string & fitsfile, bool isFront=true,
	    const std::string & extname="RPSF", size_t nrow=0);
  PsfObject(){;}
  virtual ~PsfObject() {;}
  static int findIndex(const std::vector<double> & xx, double x);
  double angularIntegral(double energy, double theta, 
			 double phi, double radius, 
			 double time) const;
  
  
private:

   static double s_energy;
   static double s_theta;
   static double s_phi;
   static double s_time;
   static const PsfObject * s_self;
   static double IPsf_coneIntegrand(double * offset);
   static void IPsf_setStaticVariables(double energy, double theta, double phi, double time, const PsfObject * self);

  // PSF parameters, energy and cos(theta) bin defs.
  std::vector<double> m_logEs;
  std::vector<double> m_energies;
  std::vector<double> m_cosths;
  std::vector<double> m_thetas;
  std::vector<std::vector<double> > m_parVectors;
  
  // PSF scaling parameters
  double m_par0;
  double m_par1;
  double m_index;
  
  // store all of the PSF parameters
  std::vector<double> m_psf_pars;

  void readFits(const std::string & fitsfile,
		const std::string & extname="RPSF",
		size_t nrow=0);
  
  void normalize_pars(double radius=90.);
  void readScaling(const std::string & fitsfile, bool isFront, 
		   const std::string & extname="PSF_SCALING_PARAMS");
  
  double evaluate(double energy, double sep, const double * pars) const;

  void getCornerPars(double energy, double theta, double & tt,
		     double & uu, std::vector<double> & cornerEnergies,
		     std::vector<size_t> & indx) const;

  double psf_base_integral(double energy, double radius, 
			   const double * pars) const;
  
  static void generateBoundaries(const std::vector<double> & x,
				 const std::vector<double> & y,
				 const std::vector<double> & values,
				 std::vector<double> & xout,
				 std::vector<double> & yout,
				 std::vector<double> & values_out, 
                                  double xlo=0, double xhi=10., 
				 double ylo=-1., double yhi=1.);
  
  double IPsf_angularIntegral(double energy, double theta, 
			      double phi, double radius, 
			      double time=0) const;

  double Psf2_psf_base_integral(double u, double gamma) const;

  double scaleFactor(double energy) const;
  double value(double separation, double energy, double theta,
	       double phi, double time) const;
  double Psf2_psf_base_function(double u, double gamma) const;
};

////////////////////////////////////////////////////////////
double PsfObject::s_energy(1e3);
double PsfObject::s_theta(0);
double PsfObject::s_phi(0);
double PsfObject::s_time(0);
const PsfObject * PsfObject::s_self(0);

PsfObject::PsfObject(std::string & fitsfile, bool isFront,
		     const std::string & extname, size_t nrow) {
  facilities::Util::expandEnvVar(&fitsfile);
  readScaling(fitsfile, isFront);
  readFits(fitsfile, extname, nrow);
  normalize_pars();
}

void PsfObject::readScaling(const std::string & fitsfile, 
			    bool isFront, 
			    const std::string & extname) {
  tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
  const tip::Table * table(fileSvc.readTable(fitsfile, extname));
  std::vector<double> values;
  st_facilities::FitsTable::getVectorData(table, "PSFSCALE", values);  
  if (isFront) {
    m_par0 = values.at(0);
    m_par1 = values.at(1);
  } else {
    m_par0 = values.at(2);
    m_par1 = values.at(3);
  }
  m_index = values.at(4);
  
  m_psf_pars.resize(values.size());
  std::copy(values.begin(), values.end(), m_psf_pars.begin());
  
  delete table;
}


void PsfObject::readFits(const std::string & fitsfile,
                     const std::string & extname, 
                     size_t nrow) {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table(fileSvc.readTable(fitsfile, extname));
   const std::vector<std::string> & validFields(table->getValidFields());

   // The first four columns *must* be "ENERG_LO", "ENERG_HI", "CTHETA_LO",
   // "CTHETA_HI", in that order.
   const char * boundsName[] = {"energ_lo", "energ_hi", 
                                "ctheta_lo", "ctheta_hi"};
   for (size_t i(0); i < 4; i++) {
      if (validFields.at(i) != boundsName[i]) {
         std::ostringstream message;
         message << "latResponse::ParTables::ParTables: "
                 << "invalid header in " << fitsfile << "  "
                 << validFields.at(i) << "  " << i;
         throw std::runtime_error(message.str());
      }
   }

   // Push boundary values onto energy and theta arrays, replicating
   // parameter values along outer boundary.

   std::vector<double> elo, ehi;
   st_facilities::FitsTable::getVectorData(table, "ENERG_LO", elo, nrow);
   st_facilities::FitsTable::getVectorData(table, "ENERG_HI", ehi, nrow);
   std::vector<double> logEs;
   for (size_t k(0); k < elo.size(); k++) {
      logEs.push_back(std::log10(std::sqrt(elo[k]*ehi[k])));
   }

   std::vector<double> mulo, muhi;
   st_facilities::FitsTable::getVectorData(table, "CTHETA_LO", mulo, nrow);
   st_facilities::FitsTable::getVectorData(table, "CTHETA_HI", muhi, nrow);
   std::vector<double> cosths;
   for (size_t i(0); i < muhi.size(); i++) {
      cosths.push_back((mulo[i] + muhi[i])/2.);
   }

   size_t par_size(elo.size()*mulo.size());

   std::vector<double> values;
   for (size_t i(4); i < validFields.size(); i++) {
      const std::string & tablename(validFields[i]);
      st_facilities::FitsTable::getVectorData(table, tablename, values, nrow);
      if (values.size() != par_size) {
         std::ostringstream message;
         message << "Parameter array size does not match "
                 << "expected size based on energy and costheta arrays "
                 << "for table " << tablename
                 << " in  " << fitsfile;
         throw std::runtime_error(message.str());
      }

      std::vector<double> my_values;
      generateBoundaries(logEs, cosths, values, 
                         m_logEs, m_cosths, my_values);

      if (i == 4) {
         m_parVectors.resize(m_logEs.size()*m_cosths.size(),
                             std::vector<double>());
      }
      for (size_t j(0); j < my_values.size(); j++) {
         m_parVectors[j].push_back(my_values[j]);
      }
   }
   for (size_t k(0); k < m_logEs.size(); k++) {
      m_energies.push_back(std::pow(10., m_logEs[k]));
   }
   for (size_t j(0); j < m_cosths.size(); j++) {
      m_thetas.push_back(std::acos(m_cosths[j])*180./M_PI);
   }
   if (m_parVectors[0].size() != 6) {
      std::ostringstream message;
      message << "Number of PSF parameters in "
              << fitsfile
              << " does no match the expected number of 6.";
      throw std::runtime_error(message.str());
   }
   delete table;
}

void PsfObject::generateBoundaries(const std::vector<double> & x,
                               const std::vector<double> & y,
                               const std::vector<double> & values,
                               std::vector<double> & xout,
                               std::vector<double> & yout,
                               std::vector<double> & values_out, 
                               double xlo, double xhi, double ylo, double yhi) {
   xout.resize(x.size() + 2);
   std::copy(x.begin(), x.end(), xout.begin() + 1);
   xout.front() = xlo;
   xout.back() = xhi;

   yout.resize(y.size() + 2);
   std::copy(y.begin(), y.end(), yout.begin() + 1);
   yout.front() = ylo;
   yout.back() = yhi;

   Array array(values, x.size());
   values_out.push_back(array(0, 0));
   for (size_t i(0); i < x.size(); i++) {
      values_out.push_back(array(0, i));
   }
   values_out.push_back(array(0, x.size()-1));
   for (size_t j(0); j < y.size(); j++) {
      values_out.push_back(array(j, 0));
      for (size_t i(0); i < x.size(); i++) {
         values_out.push_back(array(j, i));
      }
      values_out.push_back(array(j, x.size()-1));
   }
   values_out.push_back(array(y.size()-1, 0));
   for (size_t i(0); i < x.size(); i++) {
      values_out.push_back(array(y.size()-1, i));
   }
   values_out.push_back(array(y.size()-1, x.size()-1));
}
                               
double PsfObject::value(double separation, double energy, double theta,
	     double phi, double time) const {
   (void)(phi);
   (void)(time);

   double tt, uu;
   std::vector<double> cornerEnergies(4);
   std::vector<size_t> indx(4);
   getCornerPars(energy, theta, tt, uu, cornerEnergies, indx);

   double sep(separation*M_PI/180.);
   std::vector<double> yvals(4);
   for (size_t i(0); i < 4; i++) {
      yvals[i] = evaluate(cornerEnergies[i], sep, &m_parVectors[indx[i]][0]);
   }

   double my_value = st_facilities::Bilinear::evaluate(tt, uu, &yvals[0]);
   return my_value;
}

void PsfObject::IPsf_setStaticVariables(double energy, double theta, double phi, double time, const PsfObject * self) {
   s_energy = energy;
   s_theta = theta;
   s_phi = phi;
   s_time = time;
   s_self = self;
}

double PsfObject::IPsf_coneIntegrand(double * offset) {
  return s_self->value(*offset, s_energy, s_theta, s_phi, s_time)
      *std::sin(*offset*M_PI/180.)*2.*M_PI*M_PI/180.;
}

double PsfObject::IPsf_angularIntegral(double energy, 
				       double theta, double phi,
				       double radius, double time) const{
   IPsf_setStaticVariables(energy, theta, phi, time, this);
   double integral;
   double err(1e-5);
   long ierr(0);
   double zero(0);
   integral = st_facilities::GaussianQuadrature::integrate(&IPsf_coneIntegrand, zero, radius,
                                               err, ierr);
   return integral;
}

void PsfObject::normalize_pars(double radius) {
   double phi(0);
   double time(0);
   size_t indx(0);
   for (size_t j(0); j < m_thetas.size(); j++) {
      for (size_t k(0); k < m_energies.size(); k++, indx++) {
         double energy(m_energies[k]);
         double norm;
         if (energy < 120.) {
	   norm = IPsf_angularIntegral(energy, m_thetas[j], phi, 
				       radius, time);
	 } else {
	   norm = psf_base_integral(energy, radius, &m_parVectors[indx][0]);
	 }
         m_parVectors[indx][0] /= norm;
      }
   }
}

void PsfObject::getCornerPars(double energy, double theta,
                          double & tt, double & uu,
                          std::vector<double> & cornerEnergies,
                          std::vector<size_t> & indx) const {
   double logE(std::log10(energy));
   double costh(std::cos(theta*M_PI/180.));
   int i(findIndex(m_logEs, logE));
   int j(findIndex(m_cosths, costh));

   tt = (logE - m_logEs[i-1])/(m_logEs[i] - m_logEs[i-1]);
   uu = (costh - m_cosths[j-1])/(m_cosths[j] - m_cosths[j-1]);
   cornerEnergies[0] = m_energies[i-1];
   cornerEnergies[1] = m_energies[i];
   cornerEnergies[2] = m_energies[i];
   cornerEnergies[3] = m_energies[i-1];

   size_t xsize(m_energies.size());
   indx[0] = xsize*(j-1) + (i-1);
   indx[1] = xsize*(j-1) + (i);
   indx[2] = xsize*(j) + (i);
   indx[3] = xsize*(j) + (i-1);
}

int PsfObject::findIndex(const std::vector<double> & xx, double x) {
   typedef std::vector<double>::const_iterator const_iterator_t;

   const_iterator_t ix(std::upper_bound(xx.begin(), xx.end(), x));
   if (ix == xx.end() && x != xx.back()) {
      std::cout << xx.front() << "  "
                << x << "  "
                << xx.back() << std::endl;
      throw std::invalid_argument("Psf3::findIndex: x out of range");
   }
   if (x == xx.back()) {
      ix = xx.end() - 1;
   } else if (x <= xx.front()) {
      ix = xx.begin() + 1;
   }
   int i(ix - xx.begin());
   return i;
}

double PsfObject::angularIntegral(double energy, double theta, 
                              double phi, double radius, double time) const {
  if (energy < 120.) {
    double value = IPsf_angularIntegral(energy, theta, phi, radius, time);
    return value;
  }

   double tt, uu;
   std::vector<double> cornerEnergies(4);
   std::vector<size_t> indx(4);
   getCornerPars(energy, theta, tt, uu, cornerEnergies, indx);
   
   std::vector<double> yvals(4);
   for (size_t i(0); i < 4; i++) {
      yvals[i] = psf_base_integral(cornerEnergies[i], radius,
                                   &m_parVectors[indx[i]][0]);
   }
   double value = st_facilities::Bilinear::evaluate(tt, uu, &yvals[0]);
   return value;
}

double PsfObject::psf_base_integral(double energy, double radius, 
                                const double * pars) const {
   double ncore(pars[0]);
   double ntail(pars[1]);
   double score(pars[2]*scaleFactor(energy));
   double stail(pars[3]*scaleFactor(energy));
   double gcore(pars[4]);
   double gtail(pars[5]);

   double sep = radius*M_PI/180.;
   double rc = sep/score;
   double uc = rc*rc/2.;

   double rt = sep/stail;
   double ut = rt*rt/2.;
   if (gcore < 0 || gtail < 0) {
      throw std::runtime_error("gamma < 0");
   }
   return (ncore*Psf2_psf_base_integral(uc, gcore)*2.*M_PI*::sqr(score) + 
           ntail*ncore*Psf2_psf_base_integral(ut, gtail)*2.*M_PI*::sqr(stail));
}

double PsfObject::Psf2_psf_base_integral(double u, double gamma) const {
   double arg(1. + u/gamma);
   if (arg < 0) {
      std::cout << "u = " << u
                << " gamma = " << gamma
                << std::endl;
      throw std::runtime_error("neg. arg to pow");
   }
   return 1. - std::pow(arg, 1. - gamma);
}

double PsfObject::scaleFactor(double energy) const {
   double tt(std::pow(energy/100., m_index));
   return std::sqrt(::sqr(m_par0*tt) + ::sqr(m_par1));
}

//////////////////////////////////////////  

class PsfValsTool : public AlgTool, virtual public IPsfTool 
{
public:
  
  PsfValsTool( const std::string& type, 
	       const std::string& name, 
	       const IInterface* parent);
  
  virtual ~PsfValsTool() { ; }
  /// @brief Initialization of the tool.
  StatusCode initialize();
  StatusCode loadPsf(const std::string psfType);
  double computePsf(const double cl_level, const double energy,
		    const double theta, const bool isFront);
  
 private:
  std::string m_psfname;
  std::string m_filename;
  PsfObject m_psf_front;
  PsfObject m_psf_back;
};

DECLARE_TOOL_FACTORY(PsfValsTool);

PsfValsTool::PsfValsTool(const std::string& type, 
			 const std::string& name, 
			 const IInterface* parent)
  : AlgTool( type, name, parent )
{    
  // Declare additional interface
  declareInterface<IPsfTool>(this); 
  return;
}

StatusCode PsfValsTool::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;
  // Instantiate the message logger.
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "PsfValTools is initializing..," << endreq;

  return sc;
}

StatusCode PsfValsTool::loadPsf(const std::string psfName) {
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  if(psfName=="P7SOURCE_V6MC"){
    std::string psfFile_front="$(ANALYSISNTUPLEDATAPATH)/psf_P7SOURCE_V6MC_front.fits";
    m_psf_front = PsfObject(psfFile_front,true);
    std::string psfFile_back="$(ANALYSISNTUPLEDATAPATH)/psf_P7SOURCE_V6MC_back.fits";
    m_psf_back = PsfObject(psfFile_back,false);
  }
  else{
    log<<MSG::FATAL<<"Could not initialize PsfValTools, psf "<<psfName<<" not supported."<<endreq;
    return StatusCode::FAILURE;
  }
  return sc;
}

//compute containment level defined by cl_level
//for an event with pars energy,theta, isFront.
//Proceed by narrowing bisection.
double PsfValsTool::computePsf(const double cl_level, 
			       const double energy,
			       const double theta, 
			       const bool isFront){

  MsgStream log(msgSvc(), name());
  
  PsfObject psf=m_psf_back;
  if(isFront) psf=m_psf_front;

  double eps=1.e-3;//tolerance
  //low estimate for the containment radius (in degrees)
  double current_low=0.;
  //high estimate for the containment radius (in degrees)
  double current_high=90;
  double current_rad=1;//degrees
  double current_val=10.;//initial value to start while loop
  std::cout<<energy<<" "<<theta<<std::endl;
  while(std::abs(current_val-cl_level)>eps){
    current_val = psf.angularIntegral(energy,theta,0.,current_rad,0.);
    std::cout<<current_rad<<" "<<current_val<<" "<<current_low<<" "<<current_high<<std::endl;
    if(current_val>cl_level){
      current_high=current_rad;
    }
    else{
      current_low=current_rad;
    }
    current_rad=0.5*(current_low+current_high);
  }
  return current_rad;
}

double PsfObject::evaluate(double energy, double sep, const double * pars) const {
   double ncore(pars[0]);
   double ntail(pars[1]);
   double score(pars[2]*scaleFactor(energy));
   double stail(pars[3]*scaleFactor(energy));
   double gcore(pars[4]);
   double gtail(pars[5]);

   double rc = sep/score;
   double uc = rc*rc/2.;

   double rt = sep/stail;
   double ut = rt*rt/2.;
   return (ncore*Psf2_psf_base_function(uc, gcore) +
           ntail*ncore*Psf2_psf_base_function(ut, gtail));
}

double PsfObject::Psf2_psf_base_function(double u, double gamma) const {
   // ugly kluge because of sloppy programming in handoff_response
   // when setting boundaries of fit parameters for the PSF.
   if (gamma == 1) {
      gamma = 1.001;
   }
   return (1. - 1./gamma)*std::pow(1. + u/gamma, -gamma);
}
