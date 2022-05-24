#ifndef ResolutionModels_H
#define ResolutionModels_H

///////////////////////////////////////////////
// Functions for MuE FastSim
//
// G.Abbiendi  4/Dec/2018 
///////////////////////////////////////////////

#include <cmath>

namespace MuE {
  
  // Models for the RMS angular spread due to
  // multiple scattering in the material and intrinsic resolution.
  //  Functions returning theta RMS (rad) with inputs
  //  thickness in radiation lengths, p: momentum in GeV
  //
  //
  // Highland formula
  //
  inline double Highland_thrms(const double & p, const double & thickness) {
    return 0.0136/p * sqrt(thickness) * (1. + 0.038 * log(thickness));
  } 
  //
  //
  // this assumes the scattering occurs in the middle of the target (factor 0.5)
  //
  inline double Target_thrms(const double & p, const double & thickness) {
    return Highland_thrms(p, 0.5*thickness);
  }
  //
  //
  // parametrization by Antonio adding intrinsic detector resolution for an
  // apparatus made of 3 modules each one with 2 CMS CBC (for x&y), 
  // with a total length of 40cm
  //
  inline double Antonio_thrms(const double & p) {
    // resolution (mrad) = sqrt( (A/p)^2 + (B/p)^2 + C^2 ) 
    double A = 1.3571738; // mrad*GeV // Target Multiple scattering on 1 cm Be (half of it)
    double B = 1.2;       // mrad*GeV // Apparatus Multiple scattering
    double C = 0.059;     // mrad // Apparatus intrinsic resolution
    return sqrt( (A/p)*(A/p) + (B/p)*(B/p) + C*C );
  }

}

#endif
