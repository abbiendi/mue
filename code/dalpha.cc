#include <iostream>
#include <cmath>
#include <cstdio>
#include "dalpha.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////
// FORTRAN Function giving the contribution of lepton or top quark to Dalpha
// mass : mass of the lepton or top (in GeV)
// q2 : invariant momentum transfer (in GeV^2), negative for space-like
// i : 1=electron; 2=muon; 3=tau; 4=TOP quark
extern "C" double summa_(double *mass, double *q2, int *i);
////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FORTRAN subroutine (Fred Jegerlehner's Dalpha_had as a function of the scale. Version of 2017-2019)
// - names modified by Carlo just to avoid double counting with alphaQEDc19
// Inputs ---
//  de      : energy scale (in GeV) // spacelike region : de = -sqrt(-t)
//  dst2    : sin^2(theta_ew)
// Outputs ---
//  dder    : Dalpha_had
//  derrdersta : stat uncertainty on Dalpha_had
//  derrdersys : syst uncertainty on Dalpha_had
//  ddeg, derrdegsta, derrdegsys regard the Dg for the weak SU2 coupling:  g=g0/(1-Dg)
extern "C" void dhadr5x19cmcc_(double *de, double *dst2, 
			       double *dder, double *derrdersta, double *derrdersys, 
			       double *ddeg, double *derrdegsta, double *derrdegsys);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Leptonic Delta(alpha) from C.Carloni
//
double dalep(double & t, const double & invalfa0, const double & mass_e, const double & mass_mu) {
  const double pi = 3.1415926535897932384626433832795029;
  const double alpi = 1./invalfa0/pi;
  // masses: electron, muon, tau, top quark
  //  used in Carlo's code (for reweighting)
  double mass[4] = {mass_e, mass_mu, 1.77682, 175.6};
  double Sum = 0.;
  for (int i=0; i<4; ++i) {
    int ifla = i+1;
    Sum += summa_(&mass[i], &t, &ifla);
  }
  return alpi*Sum;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Interface to Fred Jegerlehner's Dalpha_had(q2) - NEWEST parameterisation -
// q2 = momentum transfer
//
double dahadFred_new(double q2) {
  double E = sqrt(std::abs(q2));
  double Q = q2<0. ? -E : E;
  double st2 = 0.2322;
  double der(0.), errdersta(0.), errdersys(0.), deg(0.), errdegsta(0.), errdegsys(0.);
  dhadr5x19cmcc_(&Q, &st2, &der, &errdersta, &errdersys, &deg, &errdegsta, &errdegsys);
  return der;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parametrizations of the Delta(alpha) hadronic
// iparam = 0 : pol2 
// iparam = 1 : Lepton-Like (one loop QED)
// iparam = 2 : modified Lepton-Like, (par[0])_LLmod = (par[0]/par[1])_LL 
//
double dahadPar(int iparam, double t, double par[2]) {

  double dahad(0);

  if (iparam==0) {
    dahad = par[0]*t+par[1]*t*t;
  }
  else if (iparam==2) {
    double squ = sqrt(1.-4.*par[1]/t);
    dahad = par[0]*par[1]*
      (-5./9.-4./3.*par[1]/t + (4./3.*par[1]*par[1]/t/t + par[1]/3./t -1./6.)*2./squ*log(std::abs((1.-squ)/(1.+squ))));
  }
  else if (iparam==1) {
    double squ = sqrt(1.-4.*par[1]/t);
    dahad = par[0] *
      (-5./9.-4./3.*par[1]/t + (4./3.*par[1]*par[1]/t/t + par[1]/3./t -1./6.)*2./squ*log(std::abs((1.-squ)/(1.+squ))));
  }
  else 
    {cout<<"\n ***ERROR: unknown dahad parametrization: iparam = "<< iparam <<endl;}

  return dahad;
}
////////////////////////////////////////////////////////////////////////////////////////////////
//
// Reweighting function
//
double reweight(int iparam, double par[2],
		uint npart, const std::array<double,11> & coeff,
		const double & invalfa0, const double & mass_e, const double & mass_mu,
		double & t24, double & t24k1, double & t24k2, double & t24k1k2) {
  
  double vp24 = 1 / (1 - dalep(t24, invalfa0, mass_e, mass_mu) - dahadPar(iparam, t24, par));

  double vp24k1 = 0;
  if (npart > 2) vp24k1 = 1 / (1 - dalep(t24k1, invalfa0, mass_e, mass_mu) - dahadPar(iparam, t24k1, par));

  double vp24k2 = 0;
  double vp24k1k2 = 0;
  if (npart > 3) {
    vp24k2 = 1 / (1 - dalep(t24k2, invalfa0, mass_e, mass_mu) - dahadPar(iparam, t24k2, par));
    vp24k1k2 = 1 / (1 - dalep(t24k1k2, invalfa0, mass_e, mass_mu) - dahadPar(iparam, t24k1k2, par));
  }

  //  printf("%.16f %.16f  %.16f  %.16f \n",vp24, vp24k1, vp24k2, vp24k1k2);

  return   coeff[0]
         + coeff[1]  * vp24*vp24
         + coeff[2]  * vp24*vp24k1
         + coeff[3]  * vp24*vp24k2
         + coeff[4]  * vp24*vp24k1k2
         + coeff[5]  * vp24k1*vp24k1
         + coeff[6]  * vp24k1*vp24k2
         + coeff[7]  * vp24k1*vp24k1k2
         + coeff[8]  * vp24k2*vp24k2
         + coeff[9]  * vp24k2*vp24k1k2
         + coeff[10] * vp24k1k2*vp24k1k2;
}

