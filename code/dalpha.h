#ifndef dalpha_H
#define dalpha_H

#include <array>

////////////////////////////////////////////////////////////////////////////////////////////////
// Leptonic Delta(alpha) from C.Carloni
//
double dalep(double & t, const double & invalfa0, const double & mass_e, const double & mass_mu);

//double dalep(double q2);

////////////////////////////////////////////////////////////////////////////////////////////////
// Interface to Fred Jegerlehner's Dalpha_had(t)
// t = space-like momentum transfer (t<0)
//
double dahadFred_new(double q2);

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parametrizations of the Delta(alpha) hadronic
// iparam = 0 : pol2 
// iparam = 1 : Lepton-Like (one loop QED)
// iparam = 2 : modified Lepton-Like, (par[0])_LLmod = (par[0]/par[1])_LL 
//
double dahadPar(int iparam, double t, double *par);

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Reweighting function
//  uses Da_had parametrization specified by dahadPar(iparam, t, par)
//
double reweight(int iparam, double par[2],
		uint npart, const std::array<double,11> & coeff,
		const double & invalfa0, const double & mass_e, const double & mass_mu,
		double & t24, double & t24k1, double & t24k2, double & t24k1k2);

#endif
