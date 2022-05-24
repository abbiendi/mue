////////////////////////////////////////////////////////////////
//
// MuE FastSim and Analysis of MESMER MC events (in root format)
//
// input: MuE configuration file
//
// G.Abbiendi  5/Jul/2018
////////////////////////////////////////////////////////////////
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <utility>
#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "MuEtree.h" 
#include "FastSim.h"
#include "Analysis.h"
#include "Utils.h"
#include "Inputs.h"
#include "TRandom.h"

using namespace std;
using namespace MuE;
  
int main(int argc, char* argv[]) {

  if (argc != 2) {
    cerr << "Usage : "<< argv[0] << " input_MuE_cfg_file \n";
    exit(100);
  }

  // set it to true to debug
  bool _debug_ = false;
  bool _debug_inp_ = true;

  MuE::One_Input cfg;
  cfg.configure(argv[1], _debug_inp_);
  //  cfg.configure("input.cfi", _debug_inp_);

  MuE::FS_Input fsi;
  fsi.configure(cfg.fastSim_ifname, _debug_inp_);

  MuE::AN_Input an_input;
  an_input.configure(cfg.analysis_ifname, _debug_inp_);

  // MC filenames to be opened are read in from list on an input file
  ifstream input_file_list(cfg.input_dirs_file);
  string root_ifname = cfg.ifname+".root";
  string line;
  vector<string> ifnames;
  vector<string> idirnames;
  vector<pair<string,Long_t> > mc_inputs;

  cerr << "looking for the following input root files : "<<endl;
  unsigned ifile = 0;
  while(getline(input_file_list, line)) {
    idirnames.push_back(line);
    string input_rootfile = line+root_ifname;
    cerr << "\t" << "file "<< ifile+1 << ": " << input_rootfile << endl;
    ifnames.push_back(input_rootfile);
    ifile++;
  }
  Long_t nfiles = ifnames.size();

  TString workdir(cfg.output_dir);
  int iexist = gSystem->mkdir(workdir);
  if (iexist != 0) {
    cerr << "ERROR: directory "<< workdir << " already exists. Please change the output name." << endl;
    exit(200);
  }
  //  gSystem->cd(workdir);
  cerr << "Writing analysis results to output directory : "<< workdir <<endl;

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 
  //gROOT->SetBatch(kFALSE); 

  // Initialize Root application
  TRint* app = new TRint("Root Application", &argc, argv);

  MuE::Setup* setup = new MuE::Setup(); 
  MuE::Event* event = new MuE::Event();
  
  TChain chain("MuEtree"); // chain event trees
  TChain parchain("MuEsetup"); // chain parameter trees

  cerr << "chaining trees... " << endl;
  cout << "Input root files: " << endl;
  unsigned ifi = 0;
  for (const auto & ifname : ifnames) {
    ++ifi;
    cout<< "File "<< ifi << ". " << ifname << endl;
    parchain.Add(ifname.data());
    chain.Add(ifname.data());
  }
  Long64_t nhdtot = parchain.GetEntries();
  if (nhdtot != nfiles) {
    cerr << "***ERROR: number of headers = "<<nhdtot
	      <<" does not correspond to the number of chained files = "<<nfiles<<endl;
    exit(300);  
  }
  cerr << "Number of input files         = " << nfiles << endl;
  Long64_t nevtot = chain.GetEntries();
  cerr << "Total number of input events  = " << nevtot << endl;

  parchain.SetBranchAddress("MuEparams", &setup);
  chain.SetBranchAddress("MuE", &event);

  if (_debug_) cerr << "parameter chain ... " << endl; 
  MuE::MCpara pargen; // parameters as read from the first tree in the chain
  MuE::MCstat sums;    // sums over the chained trees
  Double_t WmaxFinal(0);

  for (Long_t ifil=0; ifil < nfiles; ++ifil) {
    Long64_t ientry = parchain.LoadTree(ifil);
    if (ientry < 0) {
      cerr << "***ERROR in chaining MuEsetup, ifil = "<<ifil<<", ientry = "<<ientry<<endl;
      break;
    }
    parchain.GetEntry(ifil);
    if (_debug_) parchain.Show();
    if (ifil==0) {
      cout<<"\n"<<"========================================================================"<< endl;
      pargen = setup->GetMCpara();
      cout<<"MuE generator: "<< pargen.program_version << endl;
      cout<<"Run Number: "<< pargen.SampleTag <<endl;
      string strwgt = pargen.UNWGT ? "Unweighted" : "Weighted";
      cout<< strwgt << " events generation" << endl;
      cout<<"Mode : "<< pargen.Mode <<", nphotmode : "<<pargen.nphotmode <<endl;
      cout<<"radiation from muon, electron (1:on;0:off) : "<<pargen.radmuch <<" "<<pargen.radelch <<endl;
      if (!pargen.EXTBEAM) {
	cout<<"muon beam energy        = "<< pargen.Ebeam << " GeV" << endl;
	cout<<"RMS beam energy spread  = "<< pargen.EbeamRMS << " GeV" << endl;
      } else {
	cout<<"==> muon beam profile is taken from an external input"<<endl;
      }
      cout<<"muon charge             = "<< pargen.charge_mu << endl;
      cout<<"electron mass           = "<< pargen.mass_e << " GeV" << endl;
      cout<<"muon mass               = "<< pargen.mass_mu << " GeV" << endl;
      cout<<"1/alpha                 = "<< pargen.invalfa0 << endl;
      cout<<"soft photon cutoff      = "<< pargen.eps << endl;
      cout<<"photon IR mass          = "<< pargen.phmass << " GeV" <<endl;
      cout<<"hadronic running included (0:no;1:yes) : "<< pargen.ihadVP <<endl;
      cout<<"hadronic VP parameterisation           : "<< pargen.ihadVPfl <<endl; 
      cout<<"Minimum electron energy = "<< pargen.Emin_e << " GeV" << endl;
      cout<<"Min,Max electron angle  = "<< pargen.thmin_e  <<", "<< pargen.thmax_e << " (mrad)" <<endl;
      cout<<"Min,Max muon angle      = "<< pargen.thmin_mu <<", "<< pargen.thmax_mu << " (mrad)" <<endl;
      cout<<"Threshold lepton energy = "<< pargen.Ethr << " GeV" << endl;
      cout<<"Limit lepton angle      = "<< pargen.ththr << " mrad" <<endl;
      if (pargen.i_acopl == 0) cout<<"No Acoplanarity cut" << endl;
      else cout<<"Acoplanarity cut ("<<pargen.i_acopl<<") = "<< pargen.cut_acopl << " mrad" <<endl;
      if (pargen.i_elast == 0) cout<<"No Elasticity cut" << endl;
      else cout<<"Elasticity cut ("<<pargen.i_elast<<") = "<< pargen.cut_elast << " mrad" <<endl;
      cout<<"Initial normalization Wnorm = "<< setprecision(8) << pargen.Wnorm << " ub" << endl;
      if (pargen.READ_COEF) cout<<"Reweighting coefficients will be read in." << endl;
      if (pargen.isync == 0) cout<<"random number sequence is synchronized. " <<endl;
      else cout<<"random number sequence is not synchronized. " <<endl;
      cout<<"========================================================================"<< endl;
    }
    else {
      bool checkOk = CheckParameters(pargen, setup->GetMCpara());
      if (!(checkOk)) exit(400);
    }

    mc_inputs.push_back(make_pair(idirnames.at(ifil), setup->MCsums.Nwgt));

    // summary printouts for a single file
    cout<<"\n"<<"Input file : "<< ifil+1 << ". ==>> "<< ifnames.at(ifil) <<endl;
    cout<<"requested events     = "<< pargen.Nevreq << endl;
    cout<<"initial random seeds = "<< pargen.rnd_ext << " " << pargen.rnd_int << endl;
    cout<<"initial Wmax         = "<< pargen.Wmax <<endl;
    cout<<"................................................................"<< endl;
    cout<<"N generated events        = "<< setup->MCsums.Nevgen << endl;
    cout<<"N weights                 = "<< setup->MCsums.Nwgt << endl;
    cout<<"N negative weights        = "<< setup->MCsums.Nwgt_Negative << endl;
    cout<<"N weights above Wmax      = "<< setup->MCsums.Nwgt_OverMax << endl;
    cout<<"True Max weight           = "<< setup->MCsums.WmaxTrue <<endl;
    cout<<"Cross section             = "<< setup->MCsums.Xsec << " +/- " << setup->MCsums.XsecErr << endl;
    cout<<"Cross section (negative)  = "<< setup->MCsums.Xsec_Negative << " +/- " << setup->MCsums.Xsec_Negative_Err << endl;
    cout<<"Cross section (above max) = "<< setup->MCsums.Xsec_OverMax << " +/- " << setup->MCsums.Xsec_OverMax_Err << endl;

    // incremental sums
    sums.Nevgen         += setup->MCsums.Nevgen;
    sums.Nwgt           += setup->MCsums.Nwgt;
    sums.Nwgt_Negative  += setup->MCsums.Nwgt_Negative;
    sums.Swgt           += setup->MCsums.Swgt;
    sums.Swgt_Negative  += setup->MCsums.Swgt_Negative;
    sums.SQwgt          += setup->MCsums.SQwgt;
    sums.SQwgt_Negative += setup->MCsums.SQwgt_Negative;
    sums.Nwgt_OverMax   += setup->MCsums.Nwgt_OverMax ;

    // define true Wmax (over all files)
    if (setup->MCsums.WmaxTrue > WmaxFinal) WmaxFinal = setup->MCsums.WmaxTrue;
  }

  // Final quantities over all the files
  Double_t Avwgt = sums.Swgt / sums.Nwgt;
  Double_t ErrAvwgt = sqrt( (sums.SQwgt/sums.Nwgt - Avwgt*Avwgt) /sums.Nwgt );
  
  Double_t Avwgt_Negative(0); 
  Double_t ErrAvwgt_Negative(0); 
  if (sums.Nwgt_Negative > 0) {
    Avwgt_Negative = sums.Swgt_Negative / sums.Nwgt;
    ErrAvwgt_Negative = sqrt( (sums.SQwgt_Negative/sums.Nwgt - Avwgt_Negative*Avwgt_Negative) /sums.Nwgt );
  }

  Double_t sigma0 = pargen.Wnorm;
  sums.Xsec = sigma0 * Avwgt;
  sums.XsecErr = sigma0 * ErrAvwgt;
  sums.Xsec_Negative = sigma0 * Avwgt_Negative;
  sums.Xsec_Negative_Err = sigma0 * ErrAvwgt_Negative;

  // FINAL printouts
  cout<< endl;
  cout<<"================================================================"<< endl;
  cout<<"==== TOTAL STATISTICS =========================================="<< endl;
  cout<<"================================================================"<< endl;
  cout<<"N generated events        = "<< sums.Nevgen << endl;
  cout<<"N weights                 = "<< sums.Nwgt << endl;
  cout<<"N negative weights        = "<< sums.Nwgt_Negative << endl;
  cout<<"N weights above Wmax      = "<< sums.Nwgt_OverMax << endl;
  cout<<"True Max weight           = "<< WmaxFinal <<endl;
  cout<<"Cross section             = "<< sums.Xsec << " +/- " << sums.XsecErr << endl;
  cout<<"Cross section (negative)  = "<< sums.Xsec_Negative << " +/- " << sums.Xsec_Negative_Err << endl;
  //  cout<<"Cross section (above max) = "<< sums.Xsec_OverMax << " +/- " << sums.Xsec_OverMax_Err << endl;
  cout<<"Cross section (above max) =  *** TO BE DEFINED *** " << endl;
  cout<<"================================================================"<< endl;

  //.....................................................................
  // Initialise the muon and electron masses used in class Particle
  // consistently with those used by the generator
  Particle::set_mass_e(pargen.mass_e);
  Particle::set_mass_mu(pargen.mass_mu); 
  //.....................................................................
  
  // Fast Simulation //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MuE::FastSim fs(pargen, fsi, _debug_);  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ANALYSIS //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MuE::Analysis analyzer(cfg, pargen, fsi, an_input);
  analyzer.BeginJob();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////   EVENT  LOOP   //////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  //
  Long_t n_events = cfg.n_events;
  
  if (n_events <0) cerr << "\n" << "Kinematic distributions will be loaded from existing results." << endl;
  else {
    if (n_events == 0) {
      n_events = nevtot;
      if (sums.Nwgt != nevtot) 
	cerr << "***WARNING: total number of chained events ("<<nevtot
	     <<") is different from the expected number of generated events ("<<sums.Nwgt<<")"<<endl;
    } else {
      // reset normalisation to the number of analysed events
      sums.Nevgen = n_events;
      cerr << "Number of events to be analyzed = " << n_events << endl;
    }
  }

  if (_debug_) cerr << "event chain ... " << endl;

  // number of events with negligible weight (skipped)
  int zero_wgt_events = 0;

  for (Long_t iev=0; iev < n_events; ++iev) {
    Long64_t ientry = chain.LoadTree(iev);
    if (ientry <0) {
      cerr << "***ERROR in chaining MuE, event = "<<iev<<", ientry = "<<ientry<<endl;
      break;
    }
    if (iev % 1000000 == 0 ) cerr << "\n processing event : " << iev <<"\r"<< flush;
    
    chain.GetEntry(iev);
    if (_debug_)  chain.Show();

    Double_t evwgt = event->wgt_full;

    if (std::abs(evwgt) > 1e-17) {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fs.Process(*event);
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%
      analyzer.Analyze(*event, fs);
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    }
    else {
      zero_wgt_events++;
      fs.RandomNrSync();
    }
  }

  cout<<endl<< "End reading MuE events. Read "<< n_events << " events." <<endl;
  cout<<"number of zero-weight events = "<< zero_wgt_events <<endl;

  cout<<endl<<"last Random seed = "<<gRandom->GetSeed()<<endl;
  cout<<"last call to Rndm() = "<< setw(20) << setprecision(17) <<gRandom->Rndm() <<endl;
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  analyzer.EndJob(sums, fs, &mc_inputs);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  delete event;
  delete setup;

  cout << "Finished." << endl;

  if (!gROOT->IsBatch()) app->Run();

  return 0;
}
