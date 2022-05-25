#include <iostream>
#include <iomanip> 
#include <sstream>
#include "MuEtree.h"
#include "Utils.h"
#include "Mesmer.h"

using namespace std;
using namespace MuE;

MCpara::MCpara(ifstream & input_file, bool debug): MCpara() {
  // Read the header section
  string line, key, dump, str1, str2;
  istringstream stream;

  cout << "Start reading header section" <<endl;
  stream = input_line(input_file, debug);
  stream >> key;
  if (key != "<header>") {
    cout << "*** ERROR: unexpected format for header section." << endl;
    std::exit(100);
  }

  bool iok = true;

  // Sample Tag (run number)
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> SampleTag;
  iok  *= is_read_Ok(stream);
  stream = input_line(input_file, debug);
  
  // program version
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> str1 >> str2;
  iok *= is_read_Ok(stream);
  program_version = str1+" "+str2;
  
  // production details
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> process_ID >> dump >> dump >> running_on;
  //  iok *= is_read_Ok(stream);  // questa linea informativa e' inessenziale
  stream = input_line(input_file, debug);
  start_time = stream.str();

  // Number of requested events
  stream = input_line(input_file, debug);
  stream >> dump >> Nevreq;
  iok *= is_read_Ok(stream);

  // Unweighted (wgt=1) events ?  
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> UNWGT;
  iok *= is_read_Ok(stream);

  // generator Mode
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> Mode;
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> dump >> dump >> dump >> nphotmode;
  iok *= is_read_Ok(stream);
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> dump >> dump >> radmuch >> radelch;
  iok *= is_read_Ok(stream);

  // Initial seed for Random numbers 
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> rnd_ext >> rnd_int;
  iok *= is_read_Ok(stream);
  stream = input_line(input_file, debug);
  
  // Nominal Muon Beam energy
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);  
  stream >> Ebeam; 
  iok *= is_read_Ok(stream);

  // Beam energy Gaussian spread 
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> EbeamRMS;
  iok *= is_read_Ok(stream);

  // BEAM profile from external source?
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> str1;
  iok *= is_read_Ok(stream);
  EXTBEAM = (str1 == "no") ? false : true;

  // muon charge
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> charge_mu;
  iok *= is_read_Ok(stream);

  // muon mass
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> mass_mu;
  iok *= is_read_Ok(stream);

  // electron mass
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> mass_e;
  iok *= is_read_Ok(stream);

  // (fine structure constant)^-1  
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> invalfa0;  
  iok *= is_read_Ok(stream);

  // soft photon cutoff (technical parameter)
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> eps;
  iok *= is_read_Ok(stream);
  // photon mass (technical parameter)
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> dump >> phmass;
  iok *= is_read_Ok(stream);

  // MC generation cuts:
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  // Minimum Energy of outgoing electron
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> Emin_e;
  iok *= is_read_Ok(stream);
  // min and max electron angle
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> thmin_e;
  iok *= is_read_Ok(stream);
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> thmax_e;
  iok *= is_read_Ok(stream);
  // min and max muon angle
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> thmin_mu;
  iok *= is_read_Ok(stream);
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> thmax_mu;
  iok *= is_read_Ok(stream);
  // threshold min energy and max angle
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> Ethr;
  iok *= is_read_Ok(stream);
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> ththr;
  iok *= is_read_Ok(stream);
  // acoplanarity cut
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> i_acopl >> cut_acopl;
  iok *= is_read_Ok(stream);
  // elasticity cut
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> i_elast >> cut_elast;
  iok *= is_read_Ok(stream);
  stream = input_line(input_file, debug);

  // Cross section normalization factor 
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);  
  stream >> Wnorm; 
  iok *= is_read_Ok(stream);

  // Assumed maximum weight (for unweighted generation)
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);  
  stream >> Wmax; 
  iok *= is_read_Ok(stream);

  // Read Coefficients for reweighting ? 
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> READ_COEF;
  iok *= is_read_Ok(stream);
  
  // read other internal parameters
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> dump >> dump >> ihadVP;
  iok *= is_read_Ok(stream);

  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> dump >> dump >> ihadVPfl;
  iok *= is_read_Ok(stream);

  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> dump >> dump >> Nsearch;
  iok *= is_read_Ok(stream);

  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> dump >> dump >> dump >> dump >> dump >> Ndistr;
  iok *= is_read_Ok(stream);

  stream = input_line(input_file, debug);
  stream >> dump >> dump >> dump >> isync;
  iok *= is_read_Ok(stream);

  stream = input_line(input_file, debug);
  stream >> key;
  if (key == "</header>") {
      cout << "End reading header section." <<endl;
      // printout the input parameters
      Print();
      //
      if (!iok) {
	cout<<"*** ERROR: failed readout of header section."<<endl; 
	std::exit(150);
      }
  } else {
    cout << "*** ERROR: unexpected format for the header section." << endl;
    std::exit(200);
  }
}

void MCpara::InitandSetRunParams_mesmer(char* input_mesmer) {

  //Mesmer initialization
  int ierr = init_mesmer(input_mesmer);
  if (ierr !=0) cout<<"*** WARNING: MESMER initialisation error, ierr = "<<ierr <<endl;
  
  char mesmerversion[21];
  char hostname[21];
  char datetime[21];
  char RCorder[5];
  int areweighted;
  int extmubeam;
  int ivpwgts;

  mesmer_setup
    (&SampleTag,mesmerversion,hostname,datetime,&process_ID,&Nevreq,&areweighted,RCorder,
     &nphotmode,&radmuch,&radelch,&rnd_ext,&rnd_int,&Ebeam,&EbeamRMS,&extmubeam,&charge_mu,
     &mass_mu,&mass_e,&invalfa0,&Wnorm,&Wmax,&Emin_e,&thmin_e,&thmax_e,&thmin_mu,
     &thmax_mu,&Ethr,&ththr,&i_acopl,&cut_acopl,&i_elast,&cut_elast,&ivpwgts,&ihadVP,&ihadVPfl,&Nsearch,
     &Ndistr,&isync,&eps,&phmass); 

  program_version = mesmerversion;
  running_on = hostname;
  start_time = datetime;
  Mode = RCorder;
  UNWGT = areweighted == 1 ? false : true;
  EXTBEAM = extmubeam == 0 ? false : true;
  READ_COEF = ivpwgts == 1 ? true : false;

  // printout the input parameters
  Print();

  if (ierr != 0) {
    cout << "*** ERROR: MESMER initialisation error, ierr = "<<ierr <<endl;
    std::exit(250);
  }
} 

void MCpara::Print() const {
  cout<<"\n"<<"========================================================================"<< endl;
  cout<<"MuE generator: "<< program_version << endl;
  cout<<"Run Number: "<< SampleTag <<endl;
  cout<<"Number of events requested = "<< Nevreq <<endl;
  string strwgt = UNWGT ? "Unweighted" : "Weighted";
  cout<< strwgt << " events generation" << endl;
  cout<<"Mode : "<< Mode <<", nphotmode : "<<nphotmode <<endl;
  cout<<"radiation from muon, electron (1:on;0:off) : "<<radmuch <<" "<<radelch <<endl;
  if (!EXTBEAM) {
    cout<<"muon beam energy        = "<< Ebeam << " GeV" << endl;
    cout<<"RMS beam energy spread  = "<< EbeamRMS << " GeV" << endl;
  } else {
    cout<<"==> muon beam profile is taken from an external input"<<endl;
  }
  cout<<"muon charge             = "<< charge_mu << endl;
  cout<<"electron mass           = "<< mass_e << " GeV" << endl;
  cout<<"muon mass               = "<< mass_mu << " GeV" << endl;
  cout<<"soft photon cutoff      = "<< eps << endl;
  cout<<"photon mass             = "<< phmass << endl;
  cout<<"1/alpha                 = "<< invalfa0 << endl;
  cout<<"hadronic running included (0:no;1:yes) : "<< ihadVP <<endl;
  cout<<"hadronic VP parameterisation           : "<< ihadVPfl <<endl; 
  cout<<"Minimum electron energy = "<< Emin_e << " GeV" << endl;
  cout<<"Min,Max electron angle  = "<< thmin_e  <<", "<< thmax_e << " (mrad)" <<endl;
  cout<<"Min,Max muon angle      = "<< thmin_mu <<", "<< thmax_mu << " (mrad)" <<endl;
  cout<<"Threshold lepton energy = "<< Ethr << " GeV" << endl;
  cout<<"Limit lepton angle      = "<< ththr << " mrad" <<endl;
  if (i_acopl == 0) cout<<"No Acoplanarity cut" << endl;
  else cout<<"Acoplanarity cut ("<<i_acopl<<") = "<< cut_acopl << " mrad" <<endl;
  if (i_elast == 0) cout<<"No Elasticity cut" << endl;
  else cout<<"Elasticity cut ("<<i_elast<<") = "<< cut_elast << " mrad" <<endl;
  cout<<"initial random seeds    = "<< rnd_ext << " " << rnd_int << endl;
  if (isync == 0) cout<<"random number sequence will not be synchronized. " <<endl;
  else cout<<"random number sequence will be synchronized. " <<endl;
  cout<<"number of warmup shots = "<< Nsearch << endl;
  cout<<"initial Wmax            = "<< Wmax <<endl;
  cout<<"Initial normalization Wnorm = "<< setprecision(8) << Wnorm << " ub" << endl;
  if (READ_COEF) cout<<"Reweighting coefficients will be read in." << endl;
  cout<<"number of output distributions = "<< Ndistr << endl;
  cout<<"========================================================================"<< endl;
}

void MCpara::Dump() const {
      cout<<"\n"<<"========================================================================"<< endl;
      cout<<"SampleTag = "<< SampleTag <<endl;
      cout<<"program_version = "<< program_version <<endl;
      cout<<"process_ID = "<< process_ID <<", running_on = "<< running_on <<", start_time = "<< start_time <<endl;
      cout<<"Nevreq = "<< Nevreq <<", UNWGT = "<< UNWGT <<endl;
      cout<<"Mode = "<< Mode <<", nphotmode : "<<nphotmode << ", radmuch = "<<radmuch << ", radelch = "<<radelch <<endl;
      cout<<"rnd_ext = "<< rnd_ext << ", rnd_int = "<< rnd_int <<endl;
      cout<<"Ebeam = "<<Ebeam <<", EbeamRMS = "<< EbeamRMS << endl;
      cout<<"EXTBEAM = " << EXTBEAM << endl;
      cout<<"mass_mu = " << mass_mu << ", mass_e = " << mass_e << endl;
      cout<<"invalfa0 = " << invalfa0 << endl;
      cout<<"eps = " << eps <<endl;
      cout<<"photon mass = "<< phmass << endl;
      cout<<"Emin_e = "<< Emin_e <<endl;
      cout<<"thmin_e = " << thmin_e << ", thmax_e = " <<thmax_e <<endl;
      cout<<"thmin_mu = " << thmin_mu << ", thmax_mu = " <<thmax_mu <<endl;
      cout<<"Ethr = " << Ethr << ", ththr = "<< ththr <<endl;
      cout<<"i_acopl = "<< i_acopl << ", cut_acopl = "<<cut_acopl <<endl;
      cout<<"i_elast = "<< i_elast << ", cut_elast = "<<cut_elast <<endl;
      cout<<"Wnorm = "<< Wnorm <<", Wmax = "<< Wmax <<endl;
      cout<<"READ_COEF = "<< READ_COEF << endl;
      cout<<"ihadVP = "<< ihadVP <<", ihadVPfl = " << ihadVPfl <<endl;
      cout<<"Nsearch = "<< Nsearch <<endl;
      cout<<"Ndistr = " << Ndistr <<endl;
      cout<<"isync = " << isync <<endl;
      cout<<"========================================================================"<< endl;
}

MCstat::MCstat(ifstream & input_file, bool debug): MCstat() {
  
  // Read the footer section
  string line, key, dump, str1, str2;
  istringstream stream;
  cout <<endl<< "Start reading the footer section" <<endl;

  bool iok = true;

  // estimated cross section for the generated process
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Xsec >> dump >> XsecErr;
  iok *= is_read_Ok(stream);

  // total number of weights
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Nwgt;
  iok *= is_read_Ok(stream);

  // true maximum weight at the end 
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> WmaxTrue;
  iok *= is_read_Ok(stream);

  // number of weights greater than the assumed maximum (Wmax)
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Nwgt_OverMax;
  iok *= is_read_Ok(stream);

  // number of negative weights
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Nwgt_Negative;
  iok *= is_read_Ok(stream);

  // estimated bias on cross section due to weights greater than Wmax (in unweighted generation)
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Xsec_OverMax >> dump >> Xsec_OverMax_Err;
  iok *= is_read_Ok(stream);

  // estimated cross section contribution from negative weights
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Xsec_Negative >> dump >> Xsec_Negative_Err;
  iok *= is_read_Ok(stream);

  // final sum of weights and squared weights
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Swgt >> SQwgt;
  iok *= is_read_Ok(stream);

  // final sum of negative weights and squared weights
  stream = input_line(input_file, debug);
  stream = input_line(input_file, debug);
  stream >> Swgt_Negative >> SQwgt_Negative;
  iok *= is_read_Ok(stream);

  stream = input_line(input_file, debug);
  stream >> key;
  if (key == "</footer>") {
      cout << "End reading the footer section." <<endl;
      Print();

      if (!iok) {
	cout<<"*** ERROR: failed readout of footer section."<<endl; 
	std::exit(500);
      }

  } else {
    cout << "*** ERROR: unexpected format for the footer section." << endl;
    exit(600);
  }

}

void MCstat::SetEndofRun_mesmer() {
  finalize_mesmer(&Xsec, &XsecErr, &Nwgt, &Nevgen, &WmaxTrue, &Nwgt_OverMax, &Nwgt_Negative, &Xsec_OverMax, &Xsec_OverMax_Err, 
		  &Xsec_Negative, &Xsec_Negative_Err, &Swgt, &SQwgt, &Swgt_Negative, &SQwgt_Negative);
}

void MCstat::Print() const {
  cout<<"\n"<<"========================================================================"<< endl;
  cout<<"N generated events        = "<< Nevgen << endl;
  cout<<"N weights                 = "<< Nwgt << endl;
  cout<<"N negative weights        = "<< Nwgt_Negative << endl;
  cout<<"N weights above Wmax      = "<< Nwgt_OverMax << endl;
  cout<<"True Max weight           = "<< WmaxTrue <<endl;
  cout<<"Cross section             = "<< Xsec << " +/- " << XsecErr << endl;
  cout<<"Cross section (negative)  = "<< Xsec_Negative << " +/- " << Xsec_Negative_Err << endl;
  cout<<"Cross section (above max) = "<< Xsec_OverMax << " +/- " << Xsec_OverMax_Err << endl;
  cout<<"========================================================================"<< endl;
}

Double_t Particle::mass_e = 0;
Double_t Particle::mass_mu = 0;

Particle::Particle(std::istringstream & s) {
  s >> pdgId >> px >> py >> pz;
}

Double_t Particle::M() const {
  Double_t m = 0;
  if (std::abs(pdgId) == 13)      m = mass_mu; // mu+ / mu-
  else if (std::abs(pdgId) == 11) m = mass_e;  // e+ / e-
  return m;    
}

Double_t Particle::E() const {
  Double_t m = M();
  return sqrt(px*px + py*py + pz*pz + m*m);
}

void Particle::Print() const {
  cout<<"pdgId="<<pdgId
      <<", px="<<px
      <<", py="<<py
      <<", pz="<<pz<<endl;  
}

bool Event::Read(ifstream & input_file, bool read_coef, bool debug) {

// Read a NLO MC event from input stream

  istringstream stream = input_line(input_file, debug);
  string key;
  stream >> key;
  if (key == "<footer>") return false;
  else if (key != "<event>") {
    cout << "*** ERROR: unexpected format for input event (begin)" << endl;
    exit(300);
  }
  
  bool iok = true;

  // read run number
  stream = input_line(input_file, debug);
  stream >> RunNr;
  iok *= is_read_Ok(stream);

  // read event number
  stream = input_line(input_file, debug);
  stream >> EventNr;
  iok *= is_read_Ok(stream);
  
  // read number of final state particles
  stream = input_line(input_file, debug);
  UInt_t nfs;
  stream >> nfs;
  iok *= is_read_Ok(stream);
 
  // read MC weights (with full/no/leptonic running)
  stream = input_line(input_file, debug);
  stream  >> wgt_full >> wgt_norun >> wgt_lep;
  iok *= is_read_Ok(stream);
  
  // read LO and NLO weight (with full running)
  stream = input_line(input_file, debug);
  stream >> wgt_LO >> wgt_NLO;
  iok *= is_read_Ok(stream);

  // read the NC coefficients needed for reweighting
  if (read_coef) {
    stream = input_line(input_file, debug);
    for (auto & c : coef) stream >> c;
    iok *= is_read_Ok(stream);
  }
  
  // read incoming muon
  stream = input_line(input_file, debug);
  muin = MuE::Particle(stream);
  iok *= is_read_Ok(stream);

  // read final state particles
  fspart.clear();
  for (UInt_t i=0; i<nfs; ++i) {
    stream = input_line(input_file, debug);
    MuE::Particle part(stream);
    iok *= is_read_Ok(stream);
    fspart.push_back(part);
  }
  
  stream = input_line(input_file, debug);
  stream >> key;
  if (key != "</event>") {
    cout << "*** ERROR: unexpected format for input event (end)" << endl;
    exit(400);
  }

  if (!iok) cout << "*** ERROR: failed event readout, Run="<<RunNr<<", Event="<< EventNr <<endl;
  return iok;
}

int Event::GenerateEvent_mesmer(double* pmu) {

  int nfs;
  int mcids[20];
  double pmat[20][4];
  int ierr;

  generate_event_mesmer(pmu, &nfs, mcids, pmat, &wgt_full, &RunNr, &EventNr, 
			&wgt_norun, &wgt_lep, &wgt_LO, &wgt_NLO, coef.data(), &ierr);

  muin.pdgId = mcids[0];
  muin.px    = pmat[0][1];
  muin.py    = pmat[0][2];
  muin.pz    = pmat[0][3];

  MuE::Particle part;
  
  fspart.clear();
  for (int i = 0; i < nfs; ++i) {
    part.pdgId = mcids[i+2];
    part.px    = pmat[i+2][1];
    part.py    = pmat[i+2][2];
    part.pz    = pmat[i+2][3];
    fspart.push_back(part);
  }

  return ierr;
}

void Event::Print() const {
  cout<<"-------------------------------------------------"<<endl;
  cout<<"Run:"<<RunNr<<", Event:"<<EventNr<<endl;
  cout<<"wgt_full="<<std::defaultfloat<<wgt_full<<", wgt_norun="<<wgt_norun<<", wgt_lep="<<wgt_lep<<endl;
  cout<<"wgt_LO="<<wgt_LO<<", wgt_NLO="<<wgt_NLO<<endl;
  cout<<"coef=";
  for (auto & c : coef) cout<<c<<" ";
  cout<<endl;
  cout<<std::fixed<<setprecision(3)<<"incoming muon:"<<endl;
  muin.Print();
  cout<<fspart.size()<<" final state particles:"<<endl;
  for (auto & p : fspart) p.Print();
  cout<<"-------------------------------------------------"<<endl;
}
