#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Inputs.h"

using namespace MuE;
using namespace std;

void Inputs::configure(const string & cfi, bool debug) {
  ifstream ifile(cfi);
  string line;
  cerr << "Reading input cfi file: " << cfi << endl;
  getline(ifile, line);
  istringstream stream(line);
  if (debug) cerr << "\t" << stream.str() << endl;
  string key;
  stream >> key;

  if (key != "<cfi>") {
    cout << "*** ERROR: unexpected format while opening the cfi file." << endl;
    exit(100);
  }

  //--------------------
  readcfi(ifile, debug);
  //--------------------

  getline(ifile, line);
  stream = istringstream(line);
  if (debug) cerr << "\t" << stream.str() << endl;
  stream >> key;

  if (key == "</cfi>") {
    cout << "End reading cfi file." <<endl;
  } else {
    cout << "*** ERROR: unexpected format in the cfi file" << endl;
    exit(200);
  }

}

void Inputs::readin(ifstream & ifile, auto & par, bool debug) {
  string line, dump;
  getline(ifile, line);
  istringstream stream(line);
  if (debug) cerr << "\t" << stream.str() << endl;
  stream >> par >> dump;
}

void One_Input::readcfi(ifstream & ifile, bool debug) {
  readin(ifile, (*this).input_dirs_file, debug);
  readin(ifile, (*this).n_events, debug);
  readin(ifile, (*this).ifname, debug);
  readin(ifile, (*this).histo_ifname, debug);
  readin(ifile, (*this).output_dir, debug);
  readin(ifile, (*this).fastSim_ifname, debug);
  readin(ifile, (*this).analysis_ifname, debug);
}

void FS_Input::readcfi(ifstream & ifile, bool debug) {
  readin(ifile, (*this).model, debug);
  readin(ifile, (*this).MSopt, debug);
  readin(ifile, (*this).twosteps, debug);
  readin(ifile, (*this).thickness, debug);
  readin(ifile, (*this).resolution, debug);
  readin(ifile, (*this).zBias, debug);
  readin(ifile, (*this).zSigma, debug);
  if (zSigma != 0) {
    cout<<"\n ****************************************************************************"<<endl;
    cout  <<" ***WARNING: NON-ZERO zSigma breaks the random sequence in MuE ! be careful !"<<endl;
    cout  <<" ****************************************************************************"<<endl;
  }
}

void AN_Input::readcfi(ifstream & ifile, bool debug) {
  readin(ifile, (*this).makeTree, debug);
  readin(ifile, (*this).doTemplates, debug);
  readin(ifile, (*this).thetaMax, debug);
  readin(ifile, (*this).theMaxTemp, debug);
  readin(ifile, (*this).angleCorrPlots_finebins, debug);
  readin(ifile, (*this).iparam, debug);
  readin(ifile, (*this).iLumi, debug);
  readin(ifile, (*this).rangeSigma, debug);
  readin(ifile, (*this).divSigma, debug);
  readin(ifile, (*this).doEnergyScale, debug);
}
    
void AN_Input::hadrpars(double & Kref, double & Mref, double & KerrRef, double & MerrRef, double & dKu, double & dMu) const {

  cout<<"\n Parameterization of Dalpha_had(t) (0:pol2; 1:LL; 2:LLmod) = "<< iparam <<endl;  
  cout<<"\n Templates computed within +/- "<<rangeSigma<< " sigma from the expected reference values "<<endl;
  cout<<" pitch of the grid :"<<divSigma << " divisions per one sigma "<<endl;
 
  // Set reference expected central values and errors
  //
  if (iparam == 0) {
    // for pol2 parameterization
    Kref = -0.009057;  // this is  c1
    Mref = -0.01364;  //  this is c2
    KerrRef = 0.000047;  // this is c1 error
    MerrRef = 0.00064;  //  this is c2 error
    if (iLumi==1) {
      cout<<"\n ***ERROR: LowLumi not tested with iparam=0 !"<<endl;
      cout<<"\t pol2 parametrization needs to set correctly the c2 coefficient. Not foreseen presently. \n"<<endl;
      exit(1);
    }
  }
  else if (iparam == 1) {
    // for Lepton-Like param.
    Kref = 0.00721;
    Mref = 0.0525;      // GeV^2
    KerrRef = 0.00036;  // matched to the expected 1sigma in the ideal case
    MerrRef = 0.0028;   //  "   "
    if (iLumi==1) {
      KerrRef = 0.0025; // TO BE CHECKED: matched to the expected 1sigma in the ideal case for the TestRun lumi = 5/pb
      MerrRef = 0.;     // 1D param 
    }
  }
  else if (iparam ==2) {
    // for modified Lepton-Like param. K -> Kp
    Kref = 0.13724;
    Mref = 0.0525;      // GeV^2
    KerrRef = 0.00075;  // matched to the expected 1sigma in the ideal case
    MerrRef = 0.0028;   //  "   "
    if (iLumi==1) {
      KerrRef = 0.045;  // TO BE CHECKED: matched to the expected 1sigma in the ideal case for the PilotRun lumi = 3/pb
      MerrRef = 0.;     // 1D param 
    }
  }
  else {
    cout<<"\n ***ERROR: choose a parameterization !"<<endl;
    exit(1);
  }
  
  // define the grid given the chosen range and pitch 
  dKu = KerrRef/divSigma;  // matched to the expected 1sigma in the ideal case
  dMu = MerrRef/divSigma;   //    "

  cout<<" ==>>> Number of grid points on each axis = "<<ngrid() <<endl;
  cout<<"\n Expected central values and errors:  K = "<<Kref <<" +/- "<< KerrRef <<endl;
  cout<<  "                                      M = "<<Mref <<" +/- "<< MerrRef <<endl<<endl;
  
}

