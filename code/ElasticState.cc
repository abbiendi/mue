#include <cmath>
#include "ElasticState.h"

using namespace MuE;

ElasticState::ElasticState(double energy_mu_in, double mass_mu, double mass_e, double theta_e_out):
  E_(energy_mu_in), mm_(mass_mu), me_(mass_e), the_(theta_e_out)
{
  Ee_ = Ee(the_);
  Emu_ = Emu(Ee_);
  // hard-coded thmu for the_=0 (otherwise root gives NaN, while c++ is ok)
  thmu_ = the_>0 ? Thmu(Emu_) : 0.198806e-3;
  t_ = T(Ee_);
  x_ = X(t_);
}

ElasticState::ElasticState(double energy_mu_in, double mass_mu, double mass_e):
  E_(energy_mu_in), mm_(mass_mu), me_(mass_e),
  the_(0), thmu_(0), Ee_(0), Emu_(0), t_(0), x_(0)
{
  // this serves to use the LO kinematical formulas
}

inline double ElasticState::Ee(const double & theta_e_out) const {
  // beta = velocity of C.M. in the Lab frame (corresponding to the parameter "r" in EPJC)
  double beta = sqrt(1. - (mm_/E_)*(mm_/E_)) / (1. + me_/E_);
  double z  = beta * cos(theta_e_out*1e-3); // angle is in mrad
  return ((1.+z*z)/(1.-z*z)) * me_;
}

inline double ElasticState::Emu(const double & energy_e_out) const {
  return E_ + me_ - energy_e_out;
}

inline double ElasticState::Thmu(const double & energy_mu_out) const {
  double pbeam = sqrt(E_*E_ - mm_*mm_);
  double pmu = sqrt(energy_mu_out*energy_mu_out - mm_*mm_);
  double cthmu = (energy_mu_out*(E_+me_) - E_*me_ - mm_*mm_) / (pbeam*pmu);
  return 1e3*acos(cthmu); // angle is in mrad
}

inline double ElasticState::T(const double & energy_e_out) const {
  return -2.*me_*(energy_e_out - me_);
}

//inline double ElasticState::X(const double & t) const {
//  return (t + sqrt(t*t - 4*mm_*mm_*t))/(2*mm_*mm_);
//}

// calculate t given x
//inline double ElasticState::tx(const double & x) {
//  return -mm_*x/(1.-x)*mm_*x;
//}

// calculate Ee (outgoing electron energy) given t
//inline double ElasticState::Eet(const double & t) const {
//  return me_ - t/2./me_;
//}

//inline double ElasticState::theEe(const double & Ee) const {
//  double beta = sqrt(1. - (mm_/E_)*(mm_/E_)) / (1. + me_/E_);
//  return acos(sqrt((Ee-me_)/(Ee+me_)) /beta);
//}
