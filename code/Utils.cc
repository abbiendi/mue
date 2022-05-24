#include <string>
#include <sstream>
#include <iostream>
#include "Utils.h"

using namespace std;

namespace MuE {

istringstream input_line (ifstream & input_file, bool debug) {
  string line;
  getline(input_file, line);
  istringstream stream(line);
  if (debug) cout << "\t" << stream.str() << endl;
  return stream;
}

bool CheckParameters(const MCpara & pargen_0, const MCpara & pargen) 
{
  bool checkOk = true;
  
  if (pargen.program_version != pargen_0.program_version) {
    cerr << "***ERROR: inconsistent input files. program_version = "
	 <<pargen_0.program_version<<" (ref), "<<pargen.program_version<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.UNWGT != pargen_0.UNWGT) {
    cerr << "***ERROR: inconsistent input files. UNWGT = "
	 <<pargen_0.UNWGT<<" (ref), "<<pargen.UNWGT<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Mode != pargen_0.Mode) {
    cerr << "***ERROR: inconsistent input files. Mode = "
	 <<pargen_0.Mode<<" (ref), "<<pargen.Mode<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.nphotmode != pargen_0.nphotmode) {
    cerr << "***ERROR: inconsistent input files. nphotmode = "
	 <<pargen_0.nphotmode<<" (ref), "<<pargen.nphotmode<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.radmuch != pargen_0.radmuch) {
    cerr << "***ERROR: inconsistent input files. radmuch = "
	 <<pargen_0.radmuch<<" (ref), "<<pargen.radmuch<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.radelch != pargen_0.radelch) {
    cerr << "***ERROR: inconsistent input files. radelch = "
	 <<pargen_0.radelch<<" (ref), "<<pargen.radelch<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.EXTBEAM != pargen_0.EXTBEAM) {
    cerr << "***ERROR: inconsistent input files. EXTBEAM = "
	 <<pargen_0.EXTBEAM<<" (ref), "<<pargen.EXTBEAM<<" (new)" <<endl;
    checkOk = false;
  }
  if (!pargen.EXTBEAM) {
    if (pargen.Ebeam != pargen_0.Ebeam) {
      cerr << "***ERROR: inconsistent input files. Ebeam = "
	   <<pargen_0.Ebeam<<" (ref), "<<pargen.Ebeam<<" (new)" <<endl;
      checkOk = false;
    }
    if (pargen.EbeamRMS != pargen_0.EbeamRMS) {
      cerr << "***ERROR: inconsistent input files. EbeamRMS = "
	   <<pargen_0.EbeamRMS<<" (ref), "<<pargen.EbeamRMS<<" (new)" <<endl;
      checkOk = false;
    }
  }
  if (pargen.charge_mu != pargen_0.charge_mu) {
    cerr << "***ERROR: inconsistent input files. charge_mu = "
	 <<pargen_0.charge_mu<<" (ref), "<<pargen.charge_mu<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.mass_mu != pargen_0.mass_mu) {
    cerr << "***ERROR: inconsistent input files. mass_mu = "
	 <<pargen_0.mass_mu<<" (ref), "<<pargen.mass_mu<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.mass_e != pargen_0.mass_e) {
    cerr << "***ERROR: inconsistent input files. mass_e = "
	 <<pargen_0.mass_e<<" (ref), "<<pargen.mass_e<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.invalfa0 != pargen_0.invalfa0) {
    cerr << "***ERROR: inconsistent input files. invalfa0 = "
	 <<pargen_0.invalfa0<<" (ref), "<<pargen.invalfa0<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.ihadVP != pargen_0.ihadVP) {
    cerr << "***ERROR: inconsistent input files. ihadVP = "
	 <<pargen_0.ihadVP<<" (ref), "<<pargen.ihadVP<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.ihadVPfl != pargen_0.ihadVPfl) {
    cerr << "***ERROR: inconsistent input files. ihadVPfl = "
	 <<pargen_0.ihadVPfl<<" (ref), "<<pargen.ihadVPfl<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.eps != pargen_0.eps) {
    cerr << "***ERROR: inconsistent input files. eps: "
	 << pargen_0.eps << " (ref), "<< pargen.eps << " (new)" <<endl;
    cerr << "\t this may not a problem, eps is just a technical parameter, but it has to be small enough not to affect the physics results."<<endl;
    checkOk = false;
  }
  if (pargen.phmass != pargen_0.phmass) {
    cerr << "***ERROR: inconsistent input files. phmass = "
	 << pargen_0.phmass << " (ref), "<< pargen.phmass << " (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Emin_e != pargen_0.Emin_e) {
    cerr << "***ERROR: inconsistent input files. Emin_e = "
	 <<pargen_0.Emin_e<<" (ref), "<<pargen.Emin_e<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.thmin_e != pargen_0.thmin_e) {
    cerr << "***ERROR: inconsistent input files. thmin_e = "
	 <<pargen_0.thmin_e<<" (ref), "<<pargen.thmin_e<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.thmax_e != pargen_0.thmax_e) {
    cerr << "***ERROR: inconsistent input files. thmax_e = "
	 <<pargen_0.thmax_e<<" (ref), "<<pargen.thmax_e<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.thmin_mu != pargen_0.thmin_mu) {
    cerr << "***ERROR: inconsistent input files. thmin_mu = "
	 <<pargen_0.thmin_mu<<" (ref), "<<pargen.thmin_mu<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.thmax_mu != pargen_0.thmax_mu) {
    cerr << "***ERROR: inconsistent input files. thmax_mu = "
	 <<pargen_0.thmax_mu<<" (ref), "<<pargen.thmax_mu<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Ethr != pargen_0.Ethr) {
    cerr << "***ERROR: inconsistent input files. Ethr = "
	 <<pargen_0.Ethr<<" (ref), "<<pargen.Ethr<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.ththr != pargen_0.ththr) {
    cerr << "***ERROR: inconsistent input files. ththr = "
	 <<pargen_0.ththr<<" (ref), "<<pargen.ththr<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.i_acopl != pargen_0.i_acopl) {
    cerr << "***ERROR: inconsistent input files. i_acopl = "
	 <<pargen_0.i_acopl<<" (ref), "<<pargen.i_acopl<<" (new)" <<endl;
    checkOk = false;
  } else if (pargen.i_acopl != 0 && (pargen.cut_acopl != pargen_0.cut_acopl)) {
    cerr << "***ERROR: inconsistent input files. cut_acopl = "
	 <<pargen_0.cut_acopl<<" (ref), "<<pargen.cut_acopl<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.i_elast != pargen_0.i_elast) {
    cerr << "***ERROR: inconsistent input files. i_elast = "
	 <<pargen_0.i_elast<<" (ref), "<<pargen.i_elast<<" (new)" <<endl;
    checkOk = false;
  } else if (pargen.i_elast != 0 && (pargen.cut_elast != pargen_0.cut_elast)) {
    cerr << "***ERROR: inconsistent input files. cut_elast = "
	 <<pargen_0.cut_elast<<" (ref), "<<pargen.cut_elast<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Wnorm != pargen_0.Wnorm) {
    cerr << "***ERROR: inconsistent input files. Wnorm = "
	 <<pargen_0.Wnorm<<" (ref), "<<pargen.Wnorm<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Wmax != pargen_0.Wmax) {
    cerr << "***Warning: input files have different Wmax: "
	 << pargen_0.Wmax << " (ref), "<< pargen.Wmax << " (new)" <<endl;
    cerr << "\t this may not a problem, Wmax is just a technical parameter, but it should be greater or about the true Wmax known at the end."<<endl;
  }
  if (pargen.READ_COEF != pargen_0.READ_COEF) {
    cerr << "***ERROR: inconsistent input for reweighting coefficients = "
	 << pargen_0.READ_COEF << " (ref), "<< pargen.READ_COEF << " (new)" <<endl;
    checkOk = false;
  }
  if (pargen.isync != pargen_0.isync) {
    cerr << "*** Warning: input files have random number sequences not all synchronised:"<<endl;
    cerr << "\t this may not be a problem, make sure you know what you are doing."<<endl;
  }
  
  return checkOk;
}

}
