#ifndef Utils_H
#define Utils_H

///////////////////////////////////////////////
// MuE Utility functions
//
// G.Abbiendi  4/Dec/2018 
///////////////////////////////////////////////

#include <fstream>
#include <sstream>
#include "MuEtree.h"

namespace MuE {
  bool CheckParameters(const MCpara & pargen_0, const MCpara & pargen);
  std::istringstream input_line(std::ifstream & input_file, bool debug=false);
  inline bool is_read_Ok(std::istringstream & s) {return (!(s.bad() || s.fail()));}
}

#endif
