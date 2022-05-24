#include <fstream>
#include <iostream>
#include <iomanip> 
#include <cstdlib>
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "Histos.h"
#include "dalpha.h"

using namespace MuE;
using namespace std;

Histos::Histos(const MuE::MCpara & pargen, const AN_Input & paran, const TString & nameDir, const TString & namePrefix, const bool make_dirs, bool debug):
  Emin_e(pargen.Emin_e), charge_mu(pargen.charge_mu), mass_mu(pargen.mass_mu), mass_e(pargen.mass_e), invalfa0(pargen.invalfa0),
  UNWGT(pargen.UNWGT), Wnorm(pargen.Wnorm),
  angleCorrPlots_finebins(paran.angleCorrPlots_finebins), doTemplates(paran.doTemplates), thetaMax(paran.thetaMax), theMaxTemp(paran.theMaxTemp),
  iLumi(paran.iLumi), rangeSigma(paran.rangeSigma), divSigma(paran.divSigma), iparam(paran.iparam),
  doEnergyScale(paran.doEnergyScale),
  dirname(nameDir),prefix(namePrefix),makedirs(make_dirs)
{
  // init the tuple to compare with Carlo's distributions
  kine_distributions = {
    {"t13"    , "t13"   , &h_t13},
    {"t24"    , "t24"   , &h_t24},
    {"Ee"     , "eenlab", &h_Ee},
    {"the"    , "ethlab", &h_the},
    {"Emu"    , "menlab", &h_Emu},
    {"thmu"   , "mthlab", &h_thmu},
    {"egamma" , "genlab", &h_egamma},
    {"thgamma", "gthlab", &h_thgamma}
  };
  
  if (debug) {
    cerr<<"\n Testing "<< kine_distributions.size()<<" kinematical distributions in tuple:"<< endl;
    for (const auto &im : kine_distributions) cerr<< std::get<0>(im) << " " << std::get<1>(im) << " " << std::get<2>(im)  <<endl;
  }

  if (!angleCorrPlots_finebins) cout<<"\n 2D histos of angle_mu vs angle_e in fine bins will be omitted."<<endl;

  // store the sum of squares of weights for all the histos to be created
  TH1::SetDefaultSumw2();
  
  if (makedirs) {
    gFile->mkdir("gen");
    gFile->cd("gen");
    gDirectory->mkdir("kine");
    gDirectory->mkdir("ratios");
    
    gFile->cd("gen/kine");
    gDirectory->mkdir("NLO");
    gDirectory->mkdir("LO");
    gFile->cd("gen/kine/NLO");
    gDirectory->mkdir("full");
    gDirectory->mkdir("lep");
    gDirectory->mkdir("norun");
     
    gFile->mkdir("det");
    gFile->cd("det");
    gDirectory->mkdir("kine");
    gDirectory->mkdir("ratios");
    
    gFile->cd("det/kine");
    gDirectory->mkdir("NLO");
    gDirectory->mkdir("LO");
    gFile->cd("det/kine/NLO");
    gDirectory->mkdir("full");
    gDirectory->mkdir("lep");
    gDirectory->mkdir("norun");
    
    gFile->mkdir("resolution");
    gFile->cd("resolution");
    gDirectory->mkdir("NLO");
    gDirectory->mkdir("LO");

    if (doTemplates) {
      gFile->cd("gen");
      gDirectory->mkdir("templ");
      gFile->cd("gen/templ");
      gDirectory->mkdir("NLO");
      gDirectory->mkdir("LO");
      gFile->cd("gen/templ/NLO");
      gDirectory->mkdir("full");
      gDirectory->mkdir("lep");
      gDirectory->mkdir("had");
      gFile->cd("gen/templ/LO");
      gDirectory->mkdir("full");
      gDirectory->mkdir("lep");
      gDirectory->mkdir("had");
 
      gFile->cd("det");
      gDirectory->mkdir("templ");
      gFile->cd("det/templ");
      gDirectory->mkdir("NLO");
      gDirectory->mkdir("LO");
      gFile->cd("det/templ/NLO");
      gDirectory->mkdir("full");
      gDirectory->mkdir("lep");
      gDirectory->mkdir("had");
      gFile->cd("det/templ/LO");
      gDirectory->mkdir("full");
      gDirectory->mkdir("lep");
      gDirectory->mkdir("had");
    }

    if (doEnergyScale) {
      gFile->cd("gen");
      gDirectory->mkdir("escale");
      gFile->cd("gen/escale");
      gDirectory->mkdir("NLO");
      gDirectory->mkdir("LO");

      gFile->cd("det");
      gDirectory->mkdir("escale");
      gFile->cd("det/escale");
      gDirectory->mkdir("NLO");
      gDirectory->mkdir("LO");
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// FULL RUNNING (Ideal Variables) ////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("gen/kine/NLO/full");
  //
  // Kinematical Distributions with Full Running of alpha and GEN-level quantities
  hn_t13 = new TH1D(prefix+"hn_t13","t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hn_t24 = new TH1D(prefix+"hn_t24","t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hn_Ee = new TH1D(prefix+"hn_Ee","Electron energy; E_{e} (GeV); N", 100, Emin_e, 140.);
  // max electron theta: for Emin_e=0.2 GeV set to 71 mrad; for Emin_e=5 GeV to 14 mrad
  Double_t themax = Emin_e < 5.  ?  71.  :  14. ;
  hn_the = new TH1D(prefix+"hn_the","Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hn_Emu = new TH1D(prefix+"hn_Emu","Muon energy; E_{#mu} (GeV); N", 100, 0., 150.);
  hn_thmu = new TH1D(prefix+"hn_thmu","Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hn_x13 = new TH1D(prefix+"hn_x13","x_13 (muon leg); x_{13}; N",100, 0., 1.);
  hn_x24 = new TH1D(prefix+"hn_x24","x_24 (electron leg); x_{24}; N",100, 0., 1.);
  //  hn_deltaPhi = new TH1D(prefix+"hn_deltaPhi","delta(Phi) = phi(e) - phi(mu); #Delta#phi (rad); N",100,-TMath::Pi(),+TMath::Pi());
  hn_deltaPhi = new TH1D(prefix+"hn_deltaPhi","delta(Phi) = phi(e) - phi(mu); #Delta#phi (rad); N",100,-0.5,0.5);
  if (angleCorrPlots_finebins)
    hn_thmu_vs_the = new TH2D(prefix+"hn_thmu_vs_the","#theta_{#mu} vs #theta_{e}; #theta_{e} (mrad); #theta_{#mu} (mrad)",1000, 0., 30., 200, 0., 5.);
  hn_openingAngle = new TH1D(prefix+"hn_openingAngle","mu-e opening angle; opening angle (mrad); N", 142, 0., 71.);
  hn_tripleProduct = new TH1D(prefix+"hn_tripleProduct","triple product (acoplanarity); triple product; N", 200, -1., +1.);
  hn_egamma = new TH1D(prefix+"hn_egamma","Photon energy; E_{#gamma} (GeV); N", 100, 0.75, 150.);
  hn_thgamma = new TH1D(prefix+"hn_thgamma","Photon angle; #theta_{#gamma} (mrad); N", 100, 0., 40.); // in mrad
  hn_egamma_CoM = new TH1D(prefix+"hn_egamma_CoM","Photon angle CoM frame; E_{#gamma} (GeV); N", 100, 0., 0.2);
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// FULL RUNNING (Smeared Variables) //////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("det/kine/NLO/full");
  //
  // Kinematical Distributions with Full Running of alpha and Detector-level quantities
  hsn_t13 = new TH1D(prefix+"hsn_t13","DET t_13 = (pmu_in - pmu_out)^{2}; DET t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_t24 = new TH1D(prefix+"hsn_t24","DET t_24 = (pe_in - pe_out)^{2}; DET t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_Ee = new TH1D(prefix+"hsn_Ee","DET Electron energy; DET E_{e} (GeV); N", 100, Emin_e, 140.);
  hsn_the = new TH1D(prefix+"hsn_the","DET Electron angle; DET #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsn_Emu = new TH1D(prefix+"hsn_Emu","DET Muon energy; DET E_{#mu} (GeV); N", 100, 0., 150.);
  hsn_thmu = new TH1D(prefix+"hsn_thmu","DET Muon angle; DET #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hsn_x13 = new TH1D(prefix+"hsn_x13","DET x_13 (muon leg); DET x_{13}; N",100, 0., 1.);
  hsn_x24 = new TH1D(prefix+"hsn_x24","DET x_24 (electron leg); DET x_{24}; N",100, 0., 1.);
  //  hsn_deltaPhi = new TH1D(prefix+"hsn_deltaPhi","DET delta(Phi) = phi(e) - phi(mu); DET #Delta#phi (rad); N",100,-TMath::Pi(),+TMath::Pi());
  hsn_deltaPhi = new TH1D(prefix+"hsn_deltaPhi","DET delta(Phi) = phi(e) - phi(mu); DET #Delta#phi (rad); N",100,-0.5,0.5);
  if (angleCorrPlots_finebins)
    hsn_thmu_vs_the = new TH2D(prefix+"hsn_thmu_vs_the","DET #theta_{#mu} vs #theta_{e}; DET #theta_{e} (mrad); DET #theta_{#mu} (mrad)",1000, 0., 30., 200, 0., 5.);
  hsn_openingAngle = new TH1D(prefix+"hsn_openingAngle","DET mu-e opening angle; DET opening angle (mrad); N", 142, 0., 71.);
  hsn_tripleProduct = new TH1D(prefix+"hsn_tripleProduct","DET triple product (acoplanarity); DET triple product; N", 200, -1., +1.);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// NO RUNNING (Ideal variables) ///////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("gen/kine/NLO/norun");
  //
  // Kinematical Distributions with NO Running of alpha and GEN-level quantities
  hn_t13_norun = new TH1D(prefix+"hn_t13_norun","t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hn_t24_norun = new TH1D(prefix+"hn_t24_norun","t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hn_Ee_norun = new TH1D(prefix+"hn_Ee_norun","Electron energy; E_{e} (GeV); N", 100, Emin_e, 140.);
  hn_the_norun = new TH1D(prefix+"hn_the_norun","Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hn_Emu_norun = new TH1D(prefix+"hn_Emu_norun","Muon energy; E_{#mu} (GeV); N", 100, 0., 150.);
  hn_thmu_norun = new TH1D(prefix+"hn_thmu_norun","Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hn_x13_norun = new TH1D(prefix+"hn_x13_norun","x_13 (muon leg); x_{13}; N",100, 0., 1.);
  hn_x24_norun = new TH1D(prefix+"hn_x24_norun","x_24 (electron leg); x_{24}; N",100, 0., 1.);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// NO RUNNING (Smeared variables) ///////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("det/kine/NLO/norun");
  //
  // Kinematical Distributions with NO Running of alpha and DET-level quantities
  hsn_t13_norun = new TH1D(prefix+"hsn_t13_norun","DET t_13 = (pmu_in - pmu_out)^{2}; DET t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_t24_norun = new TH1D(prefix+"hsn_t24_norun","DET t_24 = (pe_in - pe_out)^{2}; DET t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_Ee_norun = new TH1D(prefix+"hsn_Ee_norun","DET Electron energy; DET E_{e} (GeV); N", 100, Emin_e, 140.);
  hsn_the_norun = new TH1D(prefix+"hsn_the_norun","DET Electron angle; DET #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsn_Emu_norun = new TH1D(prefix+"hsn_Emu_norun","DET Muon energy; DET E_{#mu} (GeV); N", 100, 0., 150.);
  hsn_thmu_norun = new TH1D(prefix+"hsn_thmu_norun","DET Muon angle; DET #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hsn_x13_norun = new TH1D(prefix+"hsn_x13_norun","DET x_13 (muon leg); DET x_{13}; N",100, 0., 1.);
  hsn_x24_norun = new TH1D(prefix+"hsn_x24_norun","DET x_24 (electron leg); DET x_{24}; N",100, 0., 1.);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// LEPTONIC RUNNING (Ideal variables) //////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("gen/kine/NLO/lep");
  //
  // Kinematical Distributions with only Leptonic Running of alpha and GEN-level quantities
  hn_t13_lep = new TH1D(prefix+"hn_t13_lep","t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hn_t24_lep = new TH1D(prefix+"hn_t24_lep","t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hn_Ee_lep = new TH1D(prefix+"hn_Ee_lep","Electron energy; E_{e} (GeV); N", 100, Emin_e, 140.);
  hn_the_lep = new TH1D(prefix+"hn_the_lep","Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hn_Emu_lep = new TH1D(prefix+"hn_Emu_lep","Muon energy; E_{#mu} (GeV); N", 100, 0., 150.);
  hn_thmu_lep = new TH1D(prefix+"hn_thmu_lep","Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hn_x13_lep = new TH1D(prefix+"hn_x13_lep","x_13 (muon leg); x_{13}; N",100, 0., 1.);
  hn_x24_lep = new TH1D(prefix+"hn_x24_lep","x_24 (electron leg); x_{24}; N",100, 0., 1.);
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// LEPTONIC RUNNING (Smeared variables) ///////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("det/kine/NLO/lep");
  //
  // Kinematical Distributions with only Leptonic Running of alpha and DET-level quantities
  hsn_t13_lep = new TH1D(prefix+"hsn_t13_lep","DET t_13 = (pmu_in - pmu_out)^{2}; DET t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_t24_lep = new TH1D(prefix+"hsn_t24_lep","DET t_24 = (pe_in - pe_out)^{2}; DET t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_Ee_lep = new TH1D(prefix+"hsn_Ee_lep","DET Electron energy; DET E_{e} (GeV); N", 100, Emin_e, 140.);
  hsn_the_lep = new TH1D(prefix+"hsn_the_lep","DET Electron angle; DET #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsn_Emu_lep = new TH1D(prefix+"hsn_Emu_lep","DET Muon energy; DET E_{#mu} (GeV); N", 100, 0., 150.);
  hsn_thmu_lep = new TH1D(prefix+"hsn_thmu_lep","DET Muon angle; DET #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hsn_x13_lep = new TH1D(prefix+"hsn_x13_lep","DET x_13 (muon leg); DET x_{13}; N",100, 0., 1.);
  hsn_x24_lep = new TH1D(prefix+"hsn_x24_lep","DET x_24 (electron leg); DET x_{24}; N",100, 0., 1.);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// LEADING ORDER (Ideal variables) /////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("gen/kine/LO");
  //
  // Kinematical Distributions with LO and GEN-level quantities (full running)
  hn_t13_LO = new TH1D(prefix+"hn_t13_LO","t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hn_t24_LO = new TH1D(prefix+"hn_t24_LO","t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hn_Ee_LO = new TH1D(prefix+"hn_Ee_LO","Electron energy; E_{e} (GeV); N", 100, Emin_e, 140.);
  hn_the_LO = new TH1D(prefix+"hn_the_LO","Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hn_Emu_LO = new TH1D(prefix+"hn_Emu_LO","Muon energy; E_{#mu} (GeV); N", 100, 0., 150.);
  hn_thmu_LO = new TH1D(prefix+"hn_thmu_LO","Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  if (angleCorrPlots_finebins)
    hn_thmu_vs_the_LO = new TH2D(prefix+"hn_thmu_vs_the_LO","#theta_{#mu} vs #theta_{e}; #theta_{e} (mrad); #theta_{#mu} (mrad)",1000, 0., 30., 200, 0., 5.);
  hn_x13_LO = new TH1D(prefix+"hn_x13_LO","x_13; x_{13}; N", 100, 0., 1.);
  hn_x24_LO = new TH1D(prefix+"hn_x24_LO","x_24; x_{24}; N", 100, 0., 1.);
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// LEADING ORDER (Smeared variables) ///////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("det/kine/LO");
  //
  // Kinematical Distributions with LO and DET-level quantities
  hsn_t13_LO = new TH1D(prefix+"hsn_t13_LO","DET t_13 = (pmu_in - pmu_out)^{2}; DET t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_t24_LO = new TH1D(prefix+"hsn_t24_LO","DET t_24 = (pe_in - pe_out)^{2}; DET t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsn_Ee_LO = new TH1D(prefix+"hsn_Ee_LO","DET Electron energy; DET E_{e} (GeV); N", 100, Emin_e, 140.);
  hsn_the_LO = new TH1D(prefix+"hsn_the_LO","DET Electron angle; DET #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsn_Emu_LO = new TH1D(prefix+"hsn_Emu_LO","DET Muon energy; DET E_{#mu} (GeV); N", 100, 0., 150.);
  hsn_thmu_LO = new TH1D(prefix+"hsn_thmu_LO","DET Muon angle; DET #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  if (angleCorrPlots_finebins)
    hsn_thmu_vs_the_LO = new TH2D(prefix+"hsn_thmu_vs_the_LO","DET #theta_{#mu} vs #theta_{e}; #theta_{e} (mrad); #theta_{#mu} (mrad)",1000, 0., 30., 200, 0., 5.);
  hsn_x13_LO = new TH1D(prefix+"hsn_x13_LO","DET x_13; x_{13}; N", 100, 0., 1.);
  hsn_x24_LO = new TH1D(prefix+"hsn_x24_LO","DET x_24; x_{24}; N", 100, 0., 1.);
 
  ////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// Angular Resolution plots ///////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  gFile->cd("resolution/NLO");
  h_resol_thex  = new TH1D(prefix+"h_resol_thex","Electron #Delta#theta_{X} (E=1GeV); #Delta#theta_{X} (mrad); N", 100,-5.,+5.);
  h_resol_they  = new TH1D(prefix+"h_resol_they","Electron #Delta#theta_{Y} (E=1GeV); #Delta#theta_{Y} (mrad); N", 100,-5.,+5.);
  h_resol_theSpaceXY  = new TH1D(prefix+"h_resol_theSpaceXY","Electron #Theta_{space} (E=1GeV); #Theta_{space} (mrad); N", 100, 0.,+5.);
  h_resol_theSpace  = new TH1D(prefix+"h_resol_theSpace","Electron #Theta_{space} (E=1GeV); #Theta_{space} (mrad); N", 100, 0.,+5.);
  h_resol_thex5  = new TH1D(prefix+"h_resol_thex5","Electron #Delta#theta_{X} (E=5GeV); #Delta#theta_{X} (mrad); N", 100,-5.,+5.);
  h_resol_they5  = new TH1D(prefix+"h_resol_they5","Electron #Delta#theta_{Y} (E=5GeV); #Delta#theta_{Y} (mrad); N", 100,-5.,+5.);
  h_resol_theSpace5XY  = new TH1D(prefix+"h_resol_theSpace5XY","Electron #Theta_{space} (E=5GeV); #Theta_{space} (mrad); N", 100, 0.,+5.);
  h_resol_theSpace5  = new TH1D(prefix+"h_resol_theSpace5","Electron #Theta_{space} (E=5GeV); #Theta_{space} (mrad); N", 100, 0.,+5.);
  h_resol_thmux = new TH1D(prefix+"h_resol_thmux","Muon #Delta#theta_{X}; #Delta#theta_{X} (mrad); N", 100,-1.,+1.);
  h_resol_thmuy = new TH1D(prefix+"h_resol_thmuy","Muon #Delta#theta_{Y}; #Delta#theta_{Y} (mrad); N", 100,-1.,+1.);
  h_resol_thmuSpaceXY = new TH1D(prefix+"h_resol_thmuSpaceXY","Muon #Theta; #Theta_{space} (mrad); N", 100, 0.,+1.);
  //
  h_resol_thexVsE  = new TH2D(prefix+"h_resol_thexVsE","Electron #Delta#theta_{X} vs Energy; E(GeV); #Delta#theta_{X}",150,0.,150., 500,-5.,+5.);
  h_resol_theyVsE  = new TH2D(prefix+"h_resol_theyVsE","Electron #Delta#theta_{Y} vs Energy; E(GeV); #Delta#theta_{Y}",150,0.,150., 500,-5.,+5.);
  h_resol_thmuxVsE = new TH2D(prefix+"h_resol_thmuxVsE","Muon #Delta#theta_{X} vs Energy; E(GeV); #Delta#theta_{X}",150,0.,150., 500,-5.,+5.);
  h_resol_thmuyVsE = new TH2D(prefix+"h_resol_thmuyVsE","Muon #Delta#theta_{Y} vs Energy; E(GeV); #Delta#theta_{Y}",150,0.,150., 500,-5.,+5.);
  //================================================================================================================
  //===================================== Ratios Smeared / Ideal Variables (with full running) =====================
  //================================================================================================================
  hsratio_t13_Resolution = new TH1D(prefix+"hsratio_t13_Resolution","t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_t24_Resolution = new TH1D(prefix+"hsratio_t24_Resolution","t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_the_Resolution = new TH1D(prefix+"hsratio_the_Resolution","Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsratio_thmu_Resolution = new TH1D(prefix+"hsratio_thmu_Resolution","Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hsratio_x13_Resolution = new TH1D(prefix+"hsratio_x13_Resolution","x_13; x_{13}; N", 100, 0., 1.);
  hsratio_x24_Resolution = new TH1D(prefix+"hsratio_x24_Resolution","x_24; x_{24}; N", 100, 0., 1.);
  //==================================================================================================================
  //===================================== Ratios Smeared / Ideal Variables (with LEADING ORDER ) =====================
  //==================================================================================================================
  gFile->cd("resolution/LO");
  hsratio_t13_Resolution_LO = new TH1D(prefix+"hsratio_t13_Resolution_LO","DET t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_t24_Resolution_LO = new TH1D(prefix+"hsratio_t24_Resolution_LO","DET t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_the_Resolution_LO = new TH1D(prefix+"hsratio_the_Resolution_LO","DET Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsratio_thmu_Resolution_LO = new TH1D(prefix+"hsratio_thmu_Resolution_LO","DET Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hsratio_x13_Resolution_LO = new TH1D(prefix+"hsratio_x13_Resolution_LO","DET x_13; x_{13}; N", 100, 0., 1.);
  hsratio_x24_Resolution_LO = new TH1D(prefix+"hsratio_x24_Resolution_LO","DET x_24; x_{24}; N", 100, 0., 1.);
 
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++ Ratios Run / NO Run = Full Running (Ideal variables) +++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  gFile->cd("gen/ratios");
  //
  hratio_t13_FullRunning = new TH1D(prefix+"hratio_t13_FullRunning","t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hratio_t24_FullRunning = new TH1D(prefix+"hratio_t24_FullRunning","t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hratio_Ee_FullRunning = new TH1D(prefix+"hratio_Ee_FullRunning","Electron energy; E_{e} (GeV); N", 100, Emin_e, 140.);
  hratio_the_FullRunning = new TH1D(prefix+"hratio_the_FullRunning","Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hratio_Emu_FullRunning = new TH1D(prefix+"hratio_Emu_FullRunning","Muon energy; E_{#mu} (GeV); N", 100, 0., 150.);
  hratio_thmu_FullRunning = new TH1D(prefix+"hratio_thmu_FullRunning","Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hratio_x13_FullRunning = new TH1D(prefix+"hratio_x13_FullRunning","x_13 (muon leg); x_{13}; N",100, 0., 1.);
  hratio_x24_FullRunning = new TH1D(prefix+"hratio_x24_FullRunning","x_24 (electron leg); x_{24}; N",100, 0., 1.);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++ Ratios Run / NO Run = Full Running (Smeared variables) +++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  gFile->cd("det/ratios");
  //
  hsratio_t13_FullRunning = new TH1D(prefix+"hsratio_t13_FullRunning","DET t_13 = (pmu_in - pmu_out)^{2}; DET t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_t24_FullRunning = new TH1D(prefix+"hsratio_t24_FullRunning","DET t_24 = (pe_in - pe_out)^{2}; DET t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_Ee_FullRunning = new TH1D(prefix+"hsratio_Ee_FullRunning","DET Electron energy; DET E_{e} (GeV); N", 100, Emin_e, 140.);
  hsratio_the_FullRunning = new TH1D(prefix+"hsratio_the_FullRunning","DET Electron angle; DET #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsratio_Emu_FullRunning = new TH1D(prefix+"hsratio_Emu_FullRunning","DET Muon energy; DET E_{#mu} (GeV); N", 100, 0., 150.);
  hsratio_thmu_FullRunning = new TH1D(prefix+"hsratio_thmu_FullRunning","DET Muon angle; DET #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hsratio_x13_FullRunning = new TH1D(prefix+"hsratio_x13_FullRunning","DET x_13 (muon leg); x_{13}; N",100, 0., 1.);
  hsratio_x24_FullRunning = new TH1D(prefix+"hsratio_x24_FullRunning","DET x_24 (electron leg); x_{24}; N",100, 0., 1.);
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++ Ratios Run / Leptonic Run = Hadronic Running (Ideal variables) +++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  gFile->cd("gen/ratios");
  //
  hratio_t13_HadronicRunning = new TH1D(prefix+"hratio_t13_HadronicRunning","t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hratio_t24_HadronicRunning = new TH1D(prefix+"hratio_t24_HadronicRunning","t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hratio_Ee_HadronicRunning = new TH1D(prefix+"hratio_Ee_HadronicRunning","Electron energy; E_{e} (GeV); N", 100, Emin_e, 140.);
  hratio_the_HadronicRunning = new TH1D(prefix+"hratio_the_HadronicRunning","Electron angle; #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hratio_Emu_HadronicRunning = new TH1D(prefix+"hratio_Emu_HadronicRunning","Muon energy; E_{#mu} (GeV); N", 100, 0., 150.);
  hratio_thmu_HadronicRunning = new TH1D(prefix+"hratio_thmu_HadronicRunning","Muon angle; #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hratio_x13_HadronicRunning = new TH1D(prefix+"hratio_x13_HadronicRunning","x_13 (muon leg); x_{13}; N",100, 0., 1.);
  hratio_x24_HadronicRunning = new TH1D(prefix+"hratio_x24_HadronicRunning","x_24 (electron leg); x_{24}; N",100, 0., 1.);
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++ Ratios Run / Leptonic Run = Hadronic Running (Smeared variables) ++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  gFile->cd("det/ratios");
  //
  hsratio_t13_HadronicRunning = new TH1D(prefix+"hsratio_t13_HadronicRunning","DET t_13 = (pmu_in - pmu_out)^{2}; DET t_{13} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_t24_HadronicRunning = new TH1D(prefix+"hsratio_t24_HadronicRunning","DET t_24 = (pe_in - pe_out)^{2}; DET t_{24} (GeV^{2}); N", 100, -0.143, 0.);
  hsratio_Ee_HadronicRunning = new TH1D(prefix+"hsratio_Ee_HadronicRunning","DET Electron energy; DET E_{e} (GeV); N", 100, Emin_e, 140.);
  hsratio_the_HadronicRunning = new TH1D(prefix+"hsratio_the_HadronicRunning","DET Electron angle; DET #theta_{e} (mrad); N", 100, 0., themax); // in mrad
  hsratio_Emu_HadronicRunning = new TH1D(prefix+"hsratio_Emu_HadronicRunning","DET Muon energy; DET E_{#mu} (GeV); N", 100, 0., 150.);
  hsratio_thmu_HadronicRunning = new TH1D(prefix+"hsratio_thmu_HadronicRunning","DET Muon angle; DET #theta_{#mu} (mrad); N", 100, 0., 5.); // in mrad
  hsratio_x13_HadronicRunning = new TH1D(prefix+"hsratio_x13_HadronicRunning","DET x_13 (muon leg); x_{13}; N",100, 0., 1.);
  hsratio_x24_HadronicRunning = new TH1D(prefix+"hsratio_x24_HadronicRunning","DET x_24 (electron leg); x_{24}; N",100, 0., 1.);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots for Energy scale calibration around equal angles
  //
  if (doEnergyScale) {

  const Int_t nbins_a = 50;
  const Double_t thetaMax_a = 5.;

  gFile->cd("det/escale/NLO");
  hsn_thA_ref = new TH1D(prefix+"hsn_thA_ref","DET Average angle; #theta_{ave} (mrad)",nbins_a, 0., thetaMax_a);
  hsn_thD_vs_thA_ref = new TH2D(prefix+"hsn_thD_vs_thA_ref","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, nbins_a, 0., thetaMax_a);
  hsn_thDs_vs_thA_ref = new TH2D(prefix+"hsn_thDs_vs_thA_ref","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, 2*nbins_a, -thetaMax_a, thetaMax_a);
  hsn_thA_ref_cut2 = new TH1D(prefix+"hsn_thA_ref_cut2","DET Average angle; #theta_{ave} (mrad)",nbins_a, 0., thetaMax_a);
  hsn_thD_vs_thA_ref_cut2 = new TH2D(prefix+"hsn_thD_vs_thA_ref_cut2","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, nbins_a, 0., thetaMax_a);
  hsn_thDs_vs_thA_ref_cut2 = new TH2D(prefix+"hsn_thDs_vs_thA_ref_cut2","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, 2*nbins_a, -thetaMax_a, thetaMax_a);

  gFile->cd("det/escale/LO");
  hsn_thA_LO_ref = new TH1D(prefix+"hsn_thA_LO_ref","DET Average angle; #theta_{ave} (mrad)",nbins_a, 0., thetaMax_a);
  hsn_thD_vs_thA_LO_ref = new TH2D(prefix+"hsn_thD_vs_thA_LO_ref","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, nbins_a, 0., thetaMax_a);
  hsn_thDs_vs_thA_LO_ref = new TH2D(prefix+"hsn_thDs_vs_thA_LO_ref","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, 2*nbins_a, -thetaMax_a, thetaMax_a);
  hsn_thA_LO_ref_cut2 = new TH1D(prefix+"hsn_thA_LO_ref_cut2","DET Average angle; #theta_{ave} (mrad)",nbins_a, 0., thetaMax_a);
  hsn_thD_vs_thA_LO_ref_cut2 = new TH2D(prefix+"hsn_thD_vs_thA_LO_ref_cut2","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, nbins_a, 0., thetaMax_a);
  hsn_thDs_vs_thA_LO_ref_cut2 = new TH2D(prefix+"hsn_thDs_vs_thA_LO_ref_cut2","DET Diff angle vs Average angle; #theta_{ave} (mrad); #theta_{dif} (mrad)",
				nbins_a, 0., thetaMax_a, 2*nbins_a, -thetaMax_a, thetaMax_a);

  //
  //=========================================================
  // Plots for Energy scale calibration (for chi2 scan)
  //
  Double_t thbin = 0.1;      // bin width: 0.1 mrad 
  Double_t thmax = thetaMax; // new v21: 50. ; old version v20: 32.;
  Int_t nbmax = thmax/thbin;
  Double_t thmumax = 6.;          // max thmu = 6 mrad
  Int_t nbmumax = thmumax/thbin;
  Double_t thavmax = thetaMax/2.; // new // old v20: 16.;
  Int_t nbavmax = thavmax/thbin;
  cout<<"\n theta histos for E scale calibration have bins of "<< thbin <<" mrad"<<endl;
  cout<<  " max theta (electron, Right, Left) = "<<thetaMax <<" mrad, (muon) = "<<thmumax<<" mrad, (average) = "<<thavmax<<" mrad"<<endl;
  cout<<  " N bins    (                     ) = "<<nbmax    <<"       (    ) = "<<nbmumax<<"       (       ) = "<<nbavmax<<endl;

  // NLO
  gFile->cd("det/escale/NLO");
  hesn_thmu_vs_the = new TH2D(prefix+"hesn_thmu_vs_the","DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbmax,0.,thmax,nbmumax,0.,thmumax);
  hesn_thL_vs_thR = new TH2D(prefix+"hesn_thL_vs_thR","DET L vs R angles; #theta_{R} (mrad); #theta_{L} (mrad)",nbmax,0.,thmax,nbmax,0.,thmax);
  hesn_thD_vs_thA = new TH2D(prefix+"hesn_thD_vs_thA","DET |#theta_{dif}| vs #theta_{ave}; #theta_{ave} (mrad); |#theta_{dif}| (mrad)",nbavmax,0.,thavmax,nbmax,0.,thmax);
  hesn_thD_vs_alpha = new TH2D(prefix+"hesn_thD_vs_alpha","DET |#theta_{dif}| vs Opening angle; #alpha (mrad); |#theta_{dif}| (mrad)",nbmax,0.,thmax,nbmax,0.,thmax);
  hesn_alpha_vs_thA = new TH2D(prefix+"hesn_alpha_vs_thA","DET Opening Angle vs #theta_{ave}; #theta_{ave} (mrad); #alpha (mrad)",nbavmax,0.,thavmax,nbmax,0.,thmax);
  // LO
  gFile->cd("det/escale/LO");
  hesn_thmu_vs_the_LO = new TH2D(prefix+"hesn_thmu_vs_the_LO","DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbmax,0.,thmax,nbmumax,0.,thmumax);
  hesn_thL_vs_thR_LO = new TH2D(prefix+"hesn_thL_vs_thR_LO","DET L vs R angles; #theta_{R} (mrad); #theta_{L} (mrad)",nbmax,0.,thmax,nbmax,0.,thmax);
  hesn_thD_vs_thA_LO = new TH2D(prefix+"hesn_thD_vs_thA_LO","DET |#theta_{dif}| vs #theta_{ave}; #theta_{ave} (mrad); |#theta_{dif}| (mrad)",nbavmax,0.,thavmax,nbmax,0.,thmax);
  hesn_thD_vs_alpha_LO = new TH2D(prefix+"hesn_thD_vs_alpha_LO","DET |#theta_{dif}| vs Opening angle; #alpha (mrad); |#theta_{dif}| (mrad)",nbmax,0.,thmax,nbmax,0.,thmax);
  hesn_alpha_vs_thA_LO = new TH2D(prefix+"hesn_alpha_vs_thA_LO","DET Opening Angle vs #theta_{ave}; #theta_{ave} (mrad); #alpha (mrad)",nbavmax,0.,thavmax,nbmax,0.,thmax);

  } //  if (doEnergyScale)

  /////////////////////////////////////////////////////////////////////////////////////////
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // T E M P L A T E  D i s t r i b u t i o n s
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  /////////////////////////////////////////////////////////////////////////////////////////
  if (!doTemplates) return;

  // calculate parameters for the hadronic running
  ngrid = paran.ngrid();
  paran.hadrpars(Kref, Mref, KerrRef, MerrRef, dKu, dMu);

  ////////////////////////////////////////////////////////////////////////
  // BINNING
  ////////////////////////////////////////////////////////////////////////
  //
  // small angle: 0 <theta <6mrad; large angle: 6mrad <theta< 32mrad
  const Double_t thetaMax_s = 6.;
  const Double_t thetaMax_l = theMaxTemp;
  // fine binning at small angle ; coarse at large angle
  const Int_t nbins_s = 60;
  const Int_t nbins_l = 26;
  const Double_t binw_s = thetaMax_s /nbins_s;                // bin width 0.1 mrad
  const Double_t binw_l = (thetaMax_l - thetaMax_s) /nbins_l; // bin width 1 mrad
  // 
  const Int_t nbins = nbins_s + nbins_l;
  //
  Double_t thetat[nbins+1]   = {0};
  Double_t thetam[nbins_s+1] = {0};
  //
  for (int i=0; i<=nbins; ++i) {
    if (i<=nbins_s) {
      thetat[i] = i*binw_s;
      thetam[i] = thetat[i];
    }
    else thetat[i] = thetaMax_s + (i-nbins_s)*binw_l;
  }
  
  //-------------------------------------------
  if (debug) {
    cout<<"\n histo total binning :"<<endl;
    for (int i=0; i<nbins+1; ++i) {
      cout<<"i = "<<i <<" "<<thetat[i] <<endl;
    }
    cout<<"\n histo muon binning :"<<endl;
    for (int i=0; i<nbins_s+1; ++i) {
      cout<<"i = "<<i <<" "<<thetam[i] <<endl;
    }
  }
  //-------------------------------------------

  //***************************************************************************************
  // Next-To-Leading Order
  //***************************************************************************************
  // theta_muon vs theta_electron  &&  theta_Left vs theta_Right  (NO ID)
  //
  // GEN-level...
  //  
  // full running Reference (with true hadronic dalfa - Jegerlehner)
  gFile->cd("gen/templ/NLO/full");
  hn_thmuVsthe_NLO_ref = new TH2D(prefix+"hn_thmuVsthe_NLO_ref","GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hn_thL_vs_thR_NLO_ref = new TH2D(prefix+"hn_thL_vs_thR_NLO_ref","GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);
  // only-leptonic running
  gFile->cd("gen/templ/NLO/lep");
  hn_thmuVsthe_NLO_lep = new TH2D(prefix+"hn_thmuVsthe_NLO_lep","GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hn_thL_vs_thR_NLO_lep = new TH2D(prefix+"hn_thL_vs_thR_NLO_lep","GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);
    
  // templates with parametrized hadronic running on the specified grid
  gFile->cd("gen/templ/NLO/had");

  if (iLumi==0) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      for (unsigned int iM=0; iM < ngrid; ++iM) {
	char name[100];
	sprintf(name,"_K%i_M%i",iK,iM);
	hn_thmuVsthe_NLO_templ.push_back(new TH2D(prefix+"hn_thmuVsthe_NLO_templ"+name,"GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
	hn_thL_vs_thR_NLO_templ.push_back(new TH2D(prefix+"hn_thL_vs_thR_NLO_templ"+name,"GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat));
      }
    }
  }
  else if (iLumi==1) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      char name[100];
      sprintf(name,"_K%i",iK);
      hn_thmuVsthe_NLO_templ.push_back(new TH2D(prefix+"hn_thmuVsthe_NLO_templ"+name,"GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
      hn_thL_vs_thR_NLO_templ.push_back(new TH2D(prefix+"hn_thL_vs_thR_NLO_templ"+name,"GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat));
    }
  }
  else {cout<<"\n *** Histos: undefined value of iLumi = "<<iLumi <<endl; exit(999);}

  //
  // DET-level...
  //  
  // full running Reference (with true hadronic dalfa - Jegerlehner)
  gFile->cd("det/templ/NLO/full");
  hsn_thmuVsthe_NLO_ref = new TH2D(prefix+"hsn_thmuVsthe_NLO_ref","DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hsn_thL_vs_thR_NLO_ref = new TH2D(prefix+"hsn_thL_vs_thR_NLO_ref","DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);
  // only-leptonic running
  gFile->cd("det/templ/NLO/lep");
  hsn_thmuVsthe_NLO_lep = new TH2D(prefix+"hsn_thmuVsthe_NLO_lep","DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hsn_thL_vs_thR_NLO_lep = new TH2D(prefix+"hsn_thL_vs_thR_NLO_lep","DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);

  // templates with parametrized hadronic running on the specified grid
  gFile->cd("det/templ/NLO/had");

  if (iLumi==0) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      for (unsigned int iM=0; iM < ngrid; ++iM) {
      char name[100];
      sprintf(name,"_K%i_M%i",iK,iM);
      hsn_thmuVsthe_NLO_templ.push_back(new TH2D(prefix+"hsn_thmuVsthe_NLO_templ"+name,"DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
      hsn_thL_vs_thR_NLO_templ.push_back(new TH2D(prefix+"hsn_thL_vs_thR_NLO_templ"+name,"DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat));
      }
    }
  }
  else if (iLumi==1) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      char name[100];
      sprintf(name,"_K%i",iK);
      hsn_thmuVsthe_NLO_templ.push_back(new TH2D(prefix+"hsn_thmuVsthe_NLO_templ"+name,"DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
      hsn_thL_vs_thR_NLO_templ.push_back(new TH2D(prefix+"hsn_thL_vs_thR_NLO_templ"+name,"DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat));
    }
  }

  //***************************************************************************************
  // Leading Order
  //***************************************************************************************
  // theta_muon vs theta_electron  &&  theta_Left vs theta_Right  (NO ID)
  //
  // GEN-level...
  //  
  // full running Reference (with true hadronic dalfa - Jegerlehner)
  gFile->cd("gen/templ/LO/full");
  hn_thmuVsthe_LO_ref = new TH2D(prefix+"hn_thmuVsthe_LO_ref","GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hn_thL_vs_thR_LO_ref = new TH2D(prefix+"hn_thL_vs_thR_LO_ref","GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);
  // only-leptonic running
  gFile->cd("gen/templ/LO/lep");
  hn_thmuVsthe_LO_lep = new TH2D(prefix+"hn_thmuVsthe_LO_lep","GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hn_thL_vs_thR_LO_lep = new TH2D(prefix+"hn_thL_vs_thR_LO_lep","GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);
    
  // templates with parametrized hadronic running on the specified grid
  gFile->cd("gen/templ/LO/had");

  if (iLumi == 0) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      for (unsigned int iM=0; iM < ngrid; ++iM) {
	char name[100];
	sprintf(name,"_K%i_M%i",iK,iM);
	hn_thmuVsthe_LO_templ.push_back(new TH2D(prefix+"hn_thmuVsthe_LO_templ"+name,"GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
	hn_thL_vs_thR_LO_templ.push_back(new TH2D(prefix+"hn_thL_vs_thR_LO_templ"+name,"GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat)); 
      }
    }
  }
  else if (iLumi == 1) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      char name[100];
      sprintf(name,"_K%i",iK);
      hn_thmuVsthe_LO_templ.push_back(new TH2D(prefix+"hn_thmuVsthe_LO_templ"+name,"GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
      hn_thL_vs_thR_LO_templ.push_back(new TH2D(prefix+"hn_thL_vs_thR_LO_templ"+name,"GEN L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat)); 
    }
  }
  
  //
  // DET-level...
  //  
  // full running Reference (with true hadronic dalfa - Jegerlehner)
  gFile->cd("det/templ/LO/full");
  hsn_thmuVsthe_LO_ref = new TH2D(prefix+"hsn_thmuVsthe_LO_ref","DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hsn_thL_vs_thR_LO_ref = new TH2D(prefix+"hsn_thL_vs_thR_LO_ref","DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);
  // only-leptonic running
  gFile->cd("det/templ/LO/lep");
  hsn_thmuVsthe_LO_lep = new TH2D(prefix+"hsn_thmuVsthe_LO_lep","DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam);
  hsn_thL_vs_thR_LO_lep = new TH2D(prefix+"hsn_thL_vs_thR_LO_lep","DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat);

  // templates with parametrized hadronic running on the specified grid
  gFile->cd("det/templ/LO/had");

  if (iLumi == 0) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      for (unsigned int iM=0; iM < ngrid; ++iM) {
	char name[100];
	sprintf(name,"_K%i_M%i",iK,iM);
	hsn_thmuVsthe_LO_templ.push_back(new TH2D(prefix+"hsn_thmuVsthe_LO_templ"+name,"DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
	hsn_thL_vs_thR_LO_templ.push_back(new TH2D(prefix+"hsn_thL_vs_thR_LO_templ"+name,"DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat));
      }
    }
  }
  else if (iLumi == 1) {
    for (unsigned int iK=0; iK < ngrid; ++iK) {
      char name[100];
      sprintf(name,"_K%i",iK);
      hsn_thmuVsthe_LO_templ.push_back(new TH2D(prefix+"hsn_thmuVsthe_LO_templ"+name,"DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)",nbins,thetat,nbins_s,thetam));
      hsn_thL_vs_thR_LO_templ.push_back(new TH2D(prefix+"hsn_thL_vs_thR_LO_templ"+name,"DET L vs R angles; #theta_{L} (mrad); #theta_{R} (mrad)",nbins,thetat,nbins,thetat));
    }
  }

}

void Histos::FillTemplates(const MuE::Event & event, const MuE::MuEana & a) {

  Double_t det_theta_L, det_theta_R;   
  Double_t gen_theta_L, gen_theta_R;

  // Define "RIGHT" the particle going most towards the X axis, i.e. with larger cos(phi), after detector simulation
  // DET variables
  if (cos(a.detKin.phmu) < cos(a.detKin.phe)) {
    det_theta_L = a.detKin.thmu;
    det_theta_R = a.detKin.the;
  } else {
    det_theta_L = a.detKin.the;
    det_theta_R = a.detKin.thmu;
  }
  // GEN variables
  if (cos(a.genKin.phmu) < cos(a.genKin.phe)) {
    gen_theta_L = a.genKin.thmu;
    gen_theta_R = a.genKin.the;
  } else {
    gen_theta_L = a.genKin.the;
    gen_theta_R = a.genKin.thmu;
  }

  // number of final state particles at GEN level
  unsigned int npart = event.fspart.size();

  // SHOULD REJECT EVENTS WITH REAL LEPTON PAIRS  (flag in MuEana ? )
  // in the moment they are exiting the execution in FastSim

  // build the templates
  //
  // calculate propagators' q2
  //
  double t24 = (a.p_e_in - a.p_e_out).M2();
  double t24k1(0);
  if (npart > 2) t24k1 = (a.p_e_in - a.p_e_out - a.V_photons[0]).M2();
  double t24k2(0);
  double t24k1k2(0);
  if (npart > 3) {
    t24k2 = (a.p_e_in - a.p_e_out - a.V_photons[1]).M2();
    t24k1k2 = (a.p_e_in - a.p_e_out - a.V_photons[0] - a.V_photons[1]).M2();
  }
  //
  // Nominal luminosity
  //----------------------
  if (iLumi == 0) {
  //----------------------
  for (unsigned int iK=0; iK < ngrid; ++iK) {
    for (unsigned int iM=0; iM < ngrid; ++iM) {
      double idK = double(iK) - double(rangeSigma*divSigma); // position in the grid in units of one division 
      double idM = double(iM) - double(rangeSigma*divSigma); //
      double K = Kref + idK*dKu; // true values of grid points 
      double M = Mref + idM*dMu; //

      double par[2] = {K, M};

      double wgtRew = reweight(iparam, par, npart, event.coef, 
			       invalfa0, mass_e, mass_mu, t24, t24k1, t24k2, t24k1k2);

      // vector index
      auto idx = ngrid*iK + iM;

      // GEN level templates
      hn_thmuVsthe_NLO_templ.at(idx)->Fill(a.genKin.the, a.genKin.thmu, wgtRew);
      hn_thL_vs_thR_NLO_templ.at(idx)->Fill(gen_theta_R, gen_theta_L, wgtRew);
      //
      // DET level templates
      hsn_thmuVsthe_NLO_templ.at(idx)->Fill(a.detKin.the, a.detKin.thmu, wgtRew);
      hsn_thL_vs_thR_NLO_templ.at(idx)->Fill(det_theta_R, det_theta_L, wgtRew);

    } // loop su M
  }   // loop su K
  //-----------------------------
  }
  // LowLumi for TestRun 2021
  //-----------------------------
  else if (iLumi == 1) {
  //-----------------------------
  for (unsigned int iK=0; iK < ngrid; ++iK) {
      double idK = double(iK) - double(rangeSigma*divSigma); // position in the grid in units of one division 
      double K = Kref + idK*dKu; // true values of grid points 
      double M = Mref; //

      double par[2] = {K, M};

      double wgtRew = reweight(iparam, par, npart, event.coef, 
			       invalfa0, mass_e, mass_mu, t24, t24k1, t24k2, t24k1k2);

      // GEN level templates
      hn_thmuVsthe_NLO_templ.at(iK)->Fill(a.genKin.the, a.genKin.thmu, wgtRew);
      hn_thL_vs_thR_NLO_templ.at(iK)->Fill(gen_theta_R, gen_theta_L, wgtRew);
      // DET level templates
      hsn_thmuVsthe_NLO_templ.at(iK)->Fill(a.detKin.the, a.detKin.thmu, wgtRew);
      hsn_thL_vs_thR_NLO_templ.at(iK)->Fill(det_theta_R, det_theta_L, wgtRew);

  }   // loop su K
  //-----------------------------
  }
  //-----------------------------

  // distributions for only-leptonic running and full running with true Jegerlehner dahad
  //
  // GEN level
  hn_thmuVsthe_NLO_lep->Fill(a.genKin.the, a.genKin.thmu, a.wgt_lep); 
  hn_thmuVsthe_NLO_ref->Fill(a.genKin.the, a.genKin.thmu, a.wgt_full);    
  hn_thL_vs_thR_NLO_lep->Fill(gen_theta_R, gen_theta_L, a.wgt_lep); 
  hn_thL_vs_thR_NLO_ref->Fill(gen_theta_R, gen_theta_L, a.wgt_full);    
  // DET level
  hsn_thmuVsthe_NLO_lep->Fill(a.detKin.the, a.detKin.thmu, a.wgt_lep);
  hsn_thmuVsthe_NLO_ref->Fill(a.detKin.the, a.detKin.thmu, a.wgt_full);
  hsn_thL_vs_thR_NLO_lep->Fill(det_theta_R, det_theta_L, a.wgt_lep); 
  hsn_thL_vs_thR_NLO_ref->Fill(det_theta_R, det_theta_L, a.wgt_full);
}

void Histos::FillTemplatesLO(const MuE::MuEana & a) {

  Double_t det_theta_L, det_theta_R;   
  Double_t gen_theta_L, gen_theta_R;

  // Define "RIGHT" the particle going most towards the X axis, i.e. with larger cos(phi), after detector simulation
  // DET variables
  if (cos(a.detKin.phmu) < cos(a.detKin.phe)) {
    det_theta_L = a.detKin.thmu;
    det_theta_R = a.detKin.the;
  } else {
    det_theta_L = a.detKin.the;
    det_theta_R = a.detKin.thmu;
  }
  // GEN variables
  if (cos(a.genKin.phmu) < cos(a.genKin.phe)) {
    gen_theta_L = a.genKin.thmu;
    gen_theta_R = a.genKin.the;
  } else {
    gen_theta_L = a.genKin.the;
    gen_theta_R = a.genKin.thmu;
  }

  // Delta(alpha) lep and had for this event
  double t = a.genKin.t13;
  double dal = dalep(t,invalfa0,mass_e,mass_mu);
  double dah = dahadFred_new(t);
  double da  = dal + dah;

  // LO weight with full running
  double wgt = a.wgt_LO; 

  // build the NLO templates
  //
  // Nominal luminosity
  //----------------------
  if (iLumi == 0) {
  //----------------------
  for (unsigned int iK=0; iK < ngrid; ++iK) {
    for (unsigned int iM=0; iM < ngrid; ++iM) {
      double idK = double(iK) - double(rangeSigma*divSigma); // position in the grid in units of one division 
      double idM = double(iM) - double(rangeSigma*divSigma); //
      double K = Kref + idK*dKu; // true values of grid points 
      double M = Mref + idM*dMu; //

      double par[2] = {K, M};

      // reweighted Da_had and Da total with the given parametrization
      double dahRew = dahadPar(iparam, t, par);
      double daRew  = dal + dahRew;
      double wgtRew = wgt * (1.-da)*(1.-da) / (1.-daRew)/(1.-daRew);

      // vector index
      auto idx = ngrid*iK + iM;

      // GEN level templates
      hn_thmuVsthe_LO_templ.at(idx)->Fill(a.genKin.the, a.genKin.thmu, wgtRew);
      hn_thL_vs_thR_LO_templ.at(idx)->Fill(gen_theta_R, gen_theta_L, wgtRew);
      // DET level templates
      hsn_thmuVsthe_LO_templ.at(idx)->Fill(a.detKin.the, a.detKin.thmu, wgtRew);
      hsn_thL_vs_thR_LO_templ.at(idx)->Fill(det_theta_R, det_theta_L, wgtRew);
    } // loop su M
  }   // loop su K
  //-----------------------------
  }
  // LowLumi for TestRun 2021
  //-----------------------------
  else if (iLumi == 1) {
  //-----------------------------
  for (unsigned int iK=0; iK < ngrid; ++iK) {
      double idK = double(iK) - double(rangeSigma*divSigma); // position in the grid in units of one division 
      double K = Kref + idK*dKu; // true values of grid points 
      double M = Mref; //

      double par[2] = {K, M};

      // reweighted Da_had and Da total with the given parametrization
      double dahRew = dahadPar(iparam, t, par);
      double daRew  = dal + dahRew;
      double wgtRew = wgt * (1.-da)*(1.-da) / (1.-daRew)/(1.-daRew);

      // GEN level templates
      hn_thmuVsthe_LO_templ.at(iK)->Fill(a.genKin.the, a.genKin.thmu, wgtRew);
      hn_thL_vs_thR_LO_templ.at(iK)->Fill(gen_theta_R, gen_theta_L, wgtRew);
      // DET level templates
      hsn_thmuVsthe_LO_templ.at(iK)->Fill(a.detKin.the, a.detKin.thmu, wgtRew);
      hsn_thL_vs_thR_LO_templ.at(iK)->Fill(det_theta_R, det_theta_L, wgtRew);
  }   // loop su K
  //-----------------------------
  }
  //-----------------------------

  // weight for leptonic-only running
  double wgtRewLep = wgt * (1.-da)*(1.-da) / (1.-dal)/(1.-dal);

  // distributions for only-leptonic running and full running with true Jegerlehner dahad
  //
  // GEN level
  hn_thmuVsthe_LO_lep->Fill(a.genKin.the, a.genKin.thmu, wgtRewLep); 
  hn_thmuVsthe_LO_ref->Fill(a.genKin.the, a.genKin.thmu, wgt);    
  hn_thL_vs_thR_LO_lep->Fill(gen_theta_R, gen_theta_L, wgtRewLep); 
  hn_thL_vs_thR_LO_ref->Fill(gen_theta_R, gen_theta_L, wgt);    
  // DET level
  hsn_thmuVsthe_LO_lep->Fill(a.detKin.the, a.detKin.thmu, wgtRewLep);
  hsn_thmuVsthe_LO_ref->Fill(a.detKin.the, a.detKin.thmu, wgt);
  hsn_thL_vs_thR_LO_lep->Fill(det_theta_R, det_theta_L, wgtRewLep); 
  hsn_thL_vs_thR_LO_ref->Fill(det_theta_R, det_theta_L, wgt);    
  
}

void Histos::FillResolution(const MuE::MuEana & a, bool debug) {

  double genTheX = atan( (tan(a.genKin.the*1e-3))*cos(a.genKin.phe) )*1e3;
  double genTheY = atan( (tan(a.genKin.the*1e-3))*sin(a.genKin.phe) )*1e3;
  double detTheX = atan( (tan(a.detKin.the*1e-3))*cos(a.detKin.phe) )*1e3;
  double detTheY = atan( (tan(a.detKin.the*1e-3))*sin(a.detKin.phe) )*1e3;
  double DTheX = detTheX - genTheX;
  double DTheY = detTheY - genTheY;
  double DTheSpace = sqrt(DTheX*DTheX + DTheY*DTheY);

  double genThmuX = atan( (tan(a.genKin.thmu*1e-3))*cos(a.genKin.phmu) )*1e3;
  double genThmuY = atan( (tan(a.genKin.thmu*1e-3))*sin(a.genKin.phmu) )*1e3;
  double detThmuX = atan( (tan(a.detKin.thmu*1e-3))*cos(a.detKin.phmu) )*1e3;
  double detThmuY = atan( (tan(a.detKin.thmu*1e-3))*sin(a.detKin.phmu) )*1e3;
  double DThmuX = detThmuX - genThmuX;
  double DThmuY = detThmuY - genThmuY;
  double DThmuSpace = sqrt(DThmuX*DThmuX + DThmuY*DThmuY);

  if (debug) {
    cout<<endl<<"Electron E = "<<a.genKin.Ee<<endl;
    cout<<"\t TheX gen = "<<genTheX <<", det = "<<detTheX <<", DTheX = "<<DTheX <<endl;
    cout<<"\t TheY gen = "<<genTheY <<", det = "<<detTheY <<", DTheX = "<<DTheY <<endl;
    cout<<endl<<"Muon E = "<<a.genKin.Emu<<endl;
    cout<<"\t ThmuX gen = "<<genThmuX <<", det = "<<detThmuX <<", DThmuX = "<<DThmuX <<endl;
    cout<<"\t ThmuY gen = "<<genThmuY <<", det = "<<detThmuY <<", DThmuX = "<<DThmuY <<endl;
  }

  //::: angle in space
  double pet = a.genKin.Ee * sin(0.001*a.genKin.the);
  double pex = pet * cos(a.genKin.phe);
  double pey = pet * sin(a.genKin.phe);
  double pez = a.genKin.Ee * cos(0.001*a.genKin.the);
  //
  double pesmt = a.detKin.Ee * sin(0.001*a.detKin.the);
  double pesmx = pesmt * cos(a.detKin.phe);
  double pesmy = pesmt * sin(a.detKin.phe);
  double pesmz = a.detKin.Ee * cos(0.001*a.detKin.the);
  //
  // from scalar product
  double cosSpace = (pex*pesmx+pey*pesmy+pez*pesmz)/(a.genKin.Ee *a.detKin.Ee); 
  double angleSpace = 1e3*acos(cosSpace);
  //:::

  if (a.genKin.Ee > 0.9 && a.genKin.Ee < 1.1) {
    h_resol_thex->Fill(DTheX, a.wgt_full);
    h_resol_they->Fill(DTheY, a.wgt_full);
    h_resol_theSpaceXY->Fill(DTheSpace, a.wgt_full);
    h_resol_theSpace->Fill(angleSpace, a.wgt_full);
  }
  if (a.genKin.Ee > 4.5 && a.genKin.Ee < 5.5) {
    h_resol_thex5->Fill(DTheX, a.wgt_full);
    h_resol_they5->Fill(DTheY, a.wgt_full);
    h_resol_theSpace5XY->Fill(DTheSpace, a.wgt_full);
    h_resol_theSpace5->Fill(angleSpace, a.wgt_full);
  }
  if (a.genKin.Emu > 148.) {
    h_resol_thmux->Fill(DThmuX, a.wgt_full);
    h_resol_thmuy->Fill(DThmuY, a.wgt_full);
    h_resol_thmuSpaceXY->Fill(DThmuSpace, a.wgt_full);
  }

  h_resol_thexVsE->Fill(a.genKin.Ee, DTheX, a.wgt_full);
  h_resol_theyVsE->Fill(a.genKin.Ee, DTheY, a.wgt_full);
  h_resol_thmuxVsE->Fill(a.genKin.Emu, DThmuX, a.wgt_full);
  h_resol_thmuyVsE->Fill(a.genKin.Emu, DThmuY, a.wgt_full);
}

void Histos::Fill(const MuE::Event & event, const MuE::MuEana & a) {
  
  // Normal Running (Leptonic + Hadronic) - GEN-level histos
  if (a.wgt_full < 1e-17) cerr<<"\n"<<"*** WARNING: wgt_full = "<< a.wgt_full
  //				   << " for Run "<<a.RunNr<< " Event "<<a.EventNr<<endl;
			      << " for Event "<<a.EventNr<<endl;
  hn_t13->Fill(a.genKin.t13, a.wgt_full);
  hn_t24->Fill(a.genKin.t24, a.wgt_full);
  hn_Ee->Fill(a.genKin.Ee, a.wgt_full);
  hn_Emu->Fill(a.genKin.Emu, a.wgt_full);
  hn_the->Fill(a.genKin.the, a.wgt_full);
  hn_thmu->Fill(a.genKin.thmu, a.wgt_full);
  hn_x13->Fill(a.genKin.x13, a.wgt_full);
  hn_x24->Fill(a.genKin.x24, a.wgt_full);
  hn_deltaPhi->Fill(a.genKin.deltaPhi, a.wgt_full);
  hn_openingAngle->Fill(a.genKin.openingAngle, a.wgt_full);
  hn_tripleProduct->Fill(a.genKin.tripleProduct, a.wgt_full);    
  //
  // Normal Running (Leptonic + Hadronic) - DET-level histos
  hsn_t13->Fill(a.detKin.t13, a.wgt_full);
  hsn_t24->Fill(a.detKin.t24, a.wgt_full);
  hsn_Ee->Fill(a.detKin.Ee, a.wgt_full);
  hsn_Emu->Fill(a.detKin.Emu, a.wgt_full);
  hsn_the->Fill(a.detKin.the, a.wgt_full);
  hsn_thmu->Fill(a.detKin.thmu, a.wgt_full);
  hsn_x13->Fill(a.detKin.x13, a.wgt_full);
  hsn_x24->Fill(a.detKin.x24, a.wgt_full);
  hsn_deltaPhi->Fill(a.detKin.deltaPhi, a.wgt_full);
  hsn_openingAngle->Fill(a.detKin.openingAngle, a.wgt_full);
  hsn_tripleProduct->Fill(a.detKin.tripleProduct, a.wgt_full);
  // 
  // 2D angle correlation in fine bins may be omitted   
  if (angleCorrPlots_finebins) {
    hn_thmu_vs_the->Fill(a.genKin.the, a.genKin.thmu, a.wgt_full);
    hsn_thmu_vs_the->Fill(a.detKin.the, a.detKin.thmu, a.wgt_full);
  }
  // radiated photon histos (GEN level)
  if (a.photon.energy >0) {
    hn_egamma->Fill(a.photon.energy, a.wgt_full);
    hn_egamma_CoM->Fill(a.photon.energyCoM, a.wgt_full);
    
    if (a.photon.energy > 0.75) {
      hn_thgamma->Fill(a.photon.theta, a.wgt_full);
    }
  }
  
  // NO running - GEN-level histos
  hn_t13_norun->Fill(a.genKin.t13, a.wgt_norun);
  hn_t24_norun->Fill(a.genKin.t24, a.wgt_norun);
  hn_Ee_norun->Fill(a.genKin.Ee, a.wgt_norun);
  hn_Emu_norun->Fill(a.genKin.Emu, a.wgt_norun);
  hn_the_norun->Fill(a.genKin.the, a.wgt_norun);
  hn_thmu_norun->Fill(a.genKin.thmu, a.wgt_norun);
  hn_x13_norun->Fill(a.genKin.x13, a.wgt_norun);
  hn_x24_norun->Fill(a.genKin.x24, a.wgt_norun);
  //    
  // NO running - DET-level histos
  hsn_t13_norun->Fill(a.detKin.t13, a.wgt_norun);
  hsn_t24_norun->Fill(a.detKin.t24, a.wgt_norun);
  hsn_Ee_norun->Fill(a.detKin.Ee, a.wgt_norun);
  hsn_Emu_norun->Fill(a.detKin.Emu, a.wgt_norun);
  hsn_the_norun->Fill(a.detKin.the, a.wgt_norun);
  hsn_thmu_norun->Fill(a.detKin.thmu, a.wgt_norun);
  hsn_x13_norun->Fill(a.detKin.x13, a.wgt_norun);
  hsn_x24_norun->Fill(a.detKin.x24, a.wgt_norun);
  
  // only Leptonic running - GEN-level histos
  hn_t13_lep->Fill(a.genKin.t13, a.wgt_lep);
  hn_t24_lep->Fill(a.genKin.t24, a.wgt_lep);
  hn_Ee_lep->Fill(a.genKin.Ee, a.wgt_lep);
  hn_Emu_lep->Fill(a.genKin.Emu, a.wgt_lep);
  hn_the_lep->Fill(a.genKin.the, a.wgt_lep);
  hn_thmu_lep->Fill(a.genKin.thmu, a.wgt_lep);
  hn_x13_lep->Fill(a.genKin.x13, a.wgt_lep);
  hn_x24_lep->Fill(a.genKin.x24, a.wgt_lep);
  //    
  // only Leptonic running - DET-level histos
  hsn_t13_lep->Fill(a.detKin.t13, a.wgt_lep);
  hsn_t24_lep->Fill(a.detKin.t24, a.wgt_lep);
  hsn_Ee_lep->Fill(a.detKin.Ee, a.wgt_lep);
  hsn_Emu_lep->Fill(a.detKin.Emu, a.wgt_lep);
  hsn_the_lep->Fill(a.detKin.the, a.wgt_lep);
  hsn_thmu_lep->Fill(a.detKin.thmu, a.wgt_lep);
  hsn_x13_lep->Fill(a.detKin.x13, a.wgt_lep);
  hsn_x24_lep->Fill(a.detKin.x24, a.wgt_lep);

  // Average angle plots (for E scale calibration)
  // --------------
  Double_t det_theta_L, det_theta_R;   
  // Define "RIGHT" the particle going most towards the X axis, i.e. with larger cos(phi), after detector simulation
  // DET variables
  if (cos(a.detKin.phmu) < cos(a.detKin.phe)) {
    det_theta_L = a.detKin.thmu;
    det_theta_R = a.detKin.the;
  } else {
    det_theta_L = a.detKin.the;
    det_theta_R = a.detKin.thmu;
  }

  Double_t theta_ave_det = (a.detKin.the + a.detKin.thmu)/2.;
  Double_t theta_dif_det = a.detKin.the - a.detKin.thmu;

  if (doEnergyScale) {
    hesn_thmu_vs_the->Fill(a.detKin.the, a.detKin.thmu, a.wgt_full);
    hesn_thL_vs_thR->Fill(det_theta_R, det_theta_L, a.wgt_full);
    hesn_thD_vs_thA->Fill(theta_ave_det, std::abs(theta_dif_det), a.wgt_full);
    hesn_thD_vs_alpha->Fill(a.detKin.openingAngle, std::abs(theta_dif_det), a.wgt_full);
    hesn_alpha_vs_thA->Fill(theta_ave_det, a.detKin.openingAngle, a.wgt_full);
    
    hsn_thA_ref->Fill(theta_ave_det, a.wgt_full);
    hsn_thD_vs_thA_ref->Fill(theta_ave_det, std::abs(theta_dif_det), a.wgt_full);
    hsn_thDs_vs_thA_ref->Fill(theta_ave_det, theta_dif_det, a.wgt_full);
    
    if (a.detKin.the > 1. && a.detKin.thmu > 1.) {
      hsn_thA_ref_cut2->Fill(theta_ave_det, a.wgt_full);
      hsn_thD_vs_thA_ref_cut2->Fill(theta_ave_det, std::abs(theta_dif_det), a.wgt_full);
      hsn_thDs_vs_thA_ref_cut2->Fill(theta_ave_det, theta_dif_det, a.wgt_full);
    }
  }

  //--------------------------------------------
  // fill templates
  if (doTemplates) FillTemplates(event, a);
  //--------------------------------------------

  // LO histos filled only when there are no photons
  if (a.wgt_LO >0.) {
    // filling with LO weights here - GEN level histos
    hn_t13_LO->Fill(a.genKin.t13, a.wgt_LO);
    hn_t24_LO->Fill(a.genKin.t24, a.wgt_LO);
    hn_Ee_LO->Fill(a.genKin.Ee, a.wgt_LO);
    hn_Emu_LO->Fill(a.genKin.Emu, a.wgt_LO);
    hn_the_LO->Fill(a.genKin.the, a.wgt_LO);
    hn_thmu_LO->Fill(a.genKin.thmu, a.wgt_LO);
    hn_x13_LO->Fill(a.genKin.x13, a.wgt_LO);
    hn_x24_LO->Fill(a.genKin.x24, a.wgt_LO);
    // LO weights - DET-level histos
    hsn_t13_LO->Fill(a.detKin.t13, a.wgt_LO);
    hsn_t24_LO->Fill(a.detKin.t24, a.wgt_LO);
    hsn_Ee_LO->Fill(a.detKin.Ee, a.wgt_LO);
    hsn_Emu_LO->Fill(a.detKin.Emu, a.wgt_LO);
    hsn_the_LO->Fill(a.detKin.the, a.wgt_LO);
    hsn_thmu_LO->Fill(a.detKin.thmu, a.wgt_LO);
    hsn_x13_LO->Fill(a.detKin.x13, a.wgt_LO);
    hsn_x24_LO->Fill(a.detKin.x24, a.wgt_LO);
    // 2D angle correlation in fine bins may be omitted   
    if (angleCorrPlots_finebins) {
      hn_thmu_vs_the_LO->Fill(a.genKin.the, a.genKin.thmu, a.wgt_LO);
      hsn_thmu_vs_the_LO->Fill(a.detKin.the, a.detKin.thmu, a.wgt_LO);
    }

    if (doEnergyScale) {
      hesn_thmu_vs_the_LO->Fill(a.detKin.the, a.detKin.thmu, a.wgt_LO);
      hesn_thL_vs_thR_LO->Fill(det_theta_R, det_theta_L, a.wgt_LO);
      hesn_thD_vs_thA_LO->Fill(theta_ave_det, std::abs(theta_dif_det), a.wgt_LO);
      hesn_thD_vs_alpha_LO->Fill(a.detKin.openingAngle, std::abs(theta_dif_det), a.wgt_LO);
      hesn_alpha_vs_thA_LO->Fill(theta_ave_det, a.detKin.openingAngle, a.wgt_LO);
      
      hsn_thA_LO_ref->Fill(theta_ave_det, a.wgt_LO);
      hsn_thD_vs_thA_LO_ref->Fill(theta_ave_det, std::abs(theta_dif_det), a.wgt_LO);
      hsn_thDs_vs_thA_LO_ref->Fill(theta_ave_det, theta_dif_det, a.wgt_LO);
      
      if (a.detKin.the > 1. && a.detKin.thmu > 1.) {
	hsn_thA_LO_ref_cut2->Fill(theta_ave_det, a.wgt_LO);
	hsn_thD_vs_thA_LO_ref_cut2->Fill(theta_ave_det, std::abs(theta_dif_det), a.wgt_LO);
	hsn_thDs_vs_thA_LO_ref_cut2->Fill(theta_ave_det, theta_dif_det, a.wgt_LO);
      }
    }

    //--------------------------------------------
    // Fill the LO TEMPLATES
    if (doTemplates) FillTemplatesLO(a);
    //--------------------------------------------
  }

  // angular resolution plots
  FillResolution(a);

}

void Histos::Normalize() {

  // normalize histograms to cross section

  if (mcsums.Nevgen <= 0) {
    cout << "*** ERROR: mcsums are empty! "<<endl;
    cout << "\t normalized histograms will not be created." <<endl;
    return;
  }

  Long_t nevtot = mcsums.Nevgen;
  Double_t wnorm = UNWGT ? mcsums.Xsec : Wnorm;
  Double_t norm = wnorm / static_cast<Double_t>(nevtot);

  // ********************************************************************
  // Next-To-Leading Order
  // ********************************************************************
  //
  // ... GEN level
  //
  gFile->cd("gen");
  if (makedirs) gDirectory->mkdir("xsection");
  gDirectory->cd("xsection");
  if (makedirs) gDirectory->mkdir("NLO");
  gDirectory->cd("NLO");

  h_t13 = (TH1D*)hn_t13->Clone(prefix+"h_t13");
  h_t13->Scale(norm, "width");
  h_t13->SetTitle("t_13 = (pmu_in - pmu_out)^{2}; t_{13} (GeV^{2}); d#sigma/dt_{13} (#mub/GeV^{2})");
  
  h_t24 = (TH1D*)hn_t24->Clone(prefix+"h_t24");
  h_t24->Scale(norm, "width");
  h_t24->SetTitle("t_24 = (pe_in - pe_out)^{2}; t_{24} (GeV^{2});  d#sigma/dt_{24} (#mub/GeV^{2})");
  
  h_Ee = (TH1D*)hn_Ee->Clone(prefix+"h_Ee");
  h_Ee->Scale(norm, "width");
  h_Ee->SetTitle("Electron energy; E_{e} (GeV); d#sigma/dE_{e} (#mub/GeV)");
  
  h_Emu = (TH1D*)hn_Emu->Clone(prefix+"h_Emu");
  h_Emu->Scale(norm, "width");
  h_Emu->SetTitle("Muon energy; E_{#mu} (GeV); d#sigma/dE_{#mu} (#mub/GeV)");
  
  h_the = (TH1D*)hn_the->Clone(prefix+"h_the");
  h_the->Scale(norm, "width");
  h_the->SetTitle("Electron angle; #theta_{e} (mrad); d#sigma/d#theta_{e} (#mub/mrad)");
  
  h_thmu = (TH1D*)hn_thmu->Clone(prefix+"h_thmu");
  h_thmu->Scale(norm, "width");
  h_thmu->SetTitle("Muon angle; #theta_{#mu} (mrad); d#sigma/d#theta_{#mu} (#mub/mrad)");
  
  h_egamma = (TH1D*)hn_egamma->Clone(prefix+"h_egamma");
  h_egamma->Scale(norm, "width");
  h_egamma->SetTitle("Photon energy; E_{#gamma} (GeV); d#sigma/dE_{#gamma} (#mub/GeV)");
  
  h_thgamma = (TH1D*)hn_thgamma->Clone(prefix+"h_thgamma");
  h_thgamma->Scale(norm, "width");
  h_thgamma->SetTitle("Photon angle; #theta_{#gamma} (mrad); d#sigma/d#theta_{#gamma} (#mub/mrad)");

  if (!doTemplates) return;
  //
  // 2D Electron-Muon angles
  //
  h_thmuVsthe_NLO_ref = (TH2D*)hn_thmuVsthe_NLO_ref->Clone(prefix+"h_thmuVsthe_NLO_ref");
  h_thmuVsthe_NLO_ref->Scale(norm, "width");
  h_thmuVsthe_NLO_ref->SetTitle("GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)");
  //
  // 2D Theta_L vs Theta_R angles  (no ID)
  //
  h_thL_vs_thR_NLO_ref = (TH2D*)hn_thL_vs_thR_NLO_ref->Clone(prefix+"h_thL_vs_thR_NLO_ref");
  h_thL_vs_thR_NLO_ref->Scale(norm, "width");
  h_thL_vs_thR_NLO_ref->SetTitle("GEN Theta_L vs Theta_R angle; #theta_{R} (mrad); #theta_{L} (mrad)");

  //
  // ... DET level
  //
  gFile->cd("det");
  if (makedirs) gDirectory->mkdir("xsection");
  gDirectory->cd("xsection");
  if (makedirs) gDirectory->mkdir("NLO");
  gDirectory->cd("NLO");
  //
  hs_thL_vs_thR_NLO_ref = (TH2D*)hsn_thL_vs_thR_NLO_ref->Clone(prefix+"hs_thL_vs_thR_NLO_ref");
  hs_thL_vs_thR_NLO_ref->Scale(norm, "width");
  hs_thL_vs_thR_NLO_ref->SetTitle("DET Theta_L vs Theta_R angle; #theta_{R} (mrad); #theta_{L} (mrad)");
  //
  hs_thmuVsthe_NLO_ref = (TH2D*)hsn_thmuVsthe_NLO_ref->Clone(prefix+"hs_thmuVsthe_NLO_ref");
  hs_thmuVsthe_NLO_ref->Scale(norm, "width");
  hs_thmuVsthe_NLO_ref->SetTitle("DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)");

  // ********************************************************************
  // Leading Order
  // ********************************************************************
  //
  // ... GEN level
  //
  gFile->cd("gen/xsection");
  if (makedirs) gDirectory->mkdir("LO");
  gDirectory->cd("LO");
  //
  // 2D Electron-Muon angles
  //
  h_thmuVsthe_LO_ref = (TH2D*)hn_thmuVsthe_LO_ref->Clone(prefix+"h_thmuVsthe_LO_ref");
  h_thmuVsthe_LO_ref->Scale(norm, "width");
  h_thmuVsthe_LO_ref->SetTitle("GEN Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)");
  //
  // 2D Theta_L vs Theta_R angles  (no ID)
  //
  h_thL_vs_thR_LO_ref = (TH2D*)hn_thL_vs_thR_LO_ref->Clone(prefix+"h_thL_vs_thR_LO_ref");
  h_thL_vs_thR_LO_ref->Scale(norm, "width");
  h_thL_vs_thR_LO_ref->SetTitle("GEN Theta_L vs Theta_R angle; #theta_{R} (mrad); #theta_{L} (mrad)");
  //
  // ... DET level
  //
  gFile->cd("det/xsection");
  if (makedirs) gDirectory->mkdir("LO");
  gDirectory->cd("LO");
  //
  hs_thmuVsthe_LO_ref = (TH2D*)hsn_thmuVsthe_LO_ref->Clone(prefix+"hs_thmuVsthe_LO_ref");
  hs_thmuVsthe_LO_ref->Scale(norm, "width");
  hs_thmuVsthe_LO_ref->SetTitle("DET Electron and Muon angles; #theta_{e} (mrad); #theta_{#mu} (mrad)");
  //
  hs_thL_vs_thR_LO_ref = (TH2D*)hsn_thL_vs_thR_LO_ref->Clone(prefix+"hs_thL_vs_thR_LO_ref");
  hs_thL_vs_thR_LO_ref->Scale(norm, "width");
  hs_thL_vs_thR_LO_ref->SetTitle("DET Theta_L vs Theta_R angle; #theta_{R} (mrad); #theta_{L} (mrad)");
}


void Histos::Fit() {

  // angular resolution Fits (from NLO distributions)
  gFile->cd("resolution/NLO");

  h_resol_thex->Fit("gaus","QLW");
  h_resol_they->Fit("gaus","QLW");
  h_resol_thex5->Fit("gaus","QLW");
  h_resol_they5->Fit("gaus","QLW");
  h_resol_thmux->Fit("gaus","QLW");
  h_resol_thmuy->Fit("gaus","QLW");

  h_resol_thexVsE->FitSlicesY(0,1,150,4,"QLW");
  h_resol_theyVsE->FitSlicesY(0,1,150,4,"QLW");
  h_resol_thmuxVsE->FitSlicesY(0,1,150,4,"QLW");
  h_resol_thmuyVsE->FitSlicesY(0,1,150,4,"QLW");
}

void Histos::PlotResolutions() {

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(11);

  TString pname = dirname+"/"+prefix;

  TCanvas c_eres1("c_eres1", "Resolution - Electron 1 GeV", 700, 500);
  c_eres1.Divide(2,1);
  c_eres1.cd(1);
  h_resol_thex->Draw();
  c_eres1.cd(2);
  h_resol_they->Draw();
  c_eres1.SaveAs(pname+"resolution_e_1GeV.png");

  TCanvas c_eres5("c_eres5", "Resolution - Electron 5 GeV", 700, 500);
  c_eres5.Divide(2,1);
  c_eres5.cd(1);
  h_resol_thex5->Draw();
  c_eres5.cd(2);
  h_resol_they5->Draw();
  c_eres5.SaveAs(pname+"resolution_e_5GeV.png");

  TCanvas c_mures("c_mures", "Resolution - Muon 148-150 GeV", 700, 500);
  c_mures.Divide(2,1);
  c_mures.cd(1);
  h_resol_thmux->Draw();
  c_mures.cd(2);
  h_resol_thmuy->Draw();
  c_mures.SaveAs(pname+"resolution_mu_148-150GeV.png");

  h_resol_thexVsE_sigma = (TH1D*)gDirectory->Get("h_resol_thexVsE_2");
  h_resol_theyVsE_sigma = (TH1D*)gDirectory->Get("h_resol_theyVsE_2");
  h_resol_thexVsE_sigma->SetTitle("Electron #sigma(#Delta#theta_{X}) vs Energy; E(GeV); #sigma(#Delta#theta_{X}) (mrad)");
  h_resol_theyVsE_sigma->SetTitle("Electron #sigma(#Delta#theta_{Y}) vs Energy; E(GeV); #sigma(#Delta#theta_{Y}) (mrad)");
  
  gStyle->SetOptStat(11);
  TCanvas c_eresVsE("c_eresVsE", "Resolution vs E - electron", 700, 700);
  c_eresVsE.Divide(1,2);
  c_eresVsE.cd(1);
  h_resol_thexVsE_sigma->Draw();
  c_eresVsE.cd(2);
  h_resol_theyVsE_sigma->Draw();
  c_eresVsE.SaveAs(pname+"resolution_vs_Energy_electron.png");

  h_resol_thmuxVsE_sigma = (TH1D*)gDirectory->Get("h_resol_thmuxVsE_2");
  h_resol_thmuyVsE_sigma = (TH1D*)gDirectory->Get("h_resol_thmuyVsE_2");

  h_resol_thmuxVsE_sigma->SetTitle("Muon #sigma(#Delta#theta_{X}) vs Energy; E(GeV); #sigma(#Delta#theta_{X}) (mrad)");
  h_resol_thmuyVsE_sigma->SetTitle("Muon #sigma(#Delta#theta_{Y}) vs Energy; E(GeV); #sigma(#Delta#theta_{Y}) (mrad)");
  TCanvas c_muresVsE("c_muresVsE", "Resolution vs E - muon", 700, 700);
  c_muresVsE.Divide(1,2);
  c_muresVsE.cd(1);
  h_resol_thmuxVsE_sigma->Draw();
  c_muresVsE.cd(2);
  h_resol_thmuyVsE_sigma->Draw();
  c_muresVsE.SaveAs(pname+"resolution_vs_Energy_muon.png");
}

void Histos::hdraw(TH1D* hist) {
  hist->SetLineColor(4);
  hist->SetLineWidth(2);
  hist->SetMarkerColor(6);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.4);
  hist->Draw();
}

void Histos::Plot() {
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetOptLogy();
  gStyle->SetOptStat(kFALSE);

  TString pname = dirname+"/"+prefix;

  TCanvas c_t13 ("c_t13", "t13", 700, 500);
  c_t13.SetLogy();
  //h_t13->GetXaxis()->SetTitle("t_{13} (GeV^2)");
  hdraw(h_t13);
  c_t13.SaveAs(pname+"t13.png");
  
  TCanvas c_t24 ("c_t24", "t24", 700, 500);
  c_t24.SetLogy();
  hdraw(h_t24);
  c_t24.SaveAs(pname+"t24.png");
  
  TCanvas c_Ee ("c_Ee", "Electron energy", 700, 500);
  c_Ee.SetLogy();
  hdraw(h_Ee);
  c_Ee.SaveAs(pname+"Ee.png");
  
  TCanvas c_Emu ("c_Emu", "Muon energy", 700, 500);
  c_Emu.SetLogy();
  hdraw(h_Emu);
  c_Emu.SaveAs(pname+"Emu.png");
  
  TCanvas c_the ("c_the", "Electron angle", 700, 500);
  hdraw(h_the);
  c_the.SaveAs(pname+"the.png");
  
  TCanvas c_thmu ("c_thmu", "Muon angle", 700, 500);
  c_thmu.SetLogy();
  hdraw(h_thmu);
  c_thmu.SaveAs(pname+"thmu.png");
  
  TCanvas c_egamma ("c_egamma", "Photon energy", 700, 500);
  c_egamma.SetLogy();
  hdraw(h_egamma);
  c_egamma.SaveAs(pname+"egamma.png");
  
  TCanvas c_thgamma ("c_thgamma", "Photon angle", 700, 500);
  hdraw(h_thgamma);
  c_thgamma.SaveAs(pname+"thgamma.png");

  TCanvas c_dphi("c_dphi", "Delta(phi) - Gen level", 700, 500);
  c_dphi.SetLogy();
  hdraw(hn_deltaPhi);
  c_dphi.SaveAs(pname+"deltaPhi_gen.png");

  TCanvas cs_dphi("cs_dphi", "Delta(phi) - Det level", 700, 500);
  cs_dphi.SetLogy();
  hdraw(hsn_deltaPhi);
  cs_dphi.SaveAs(pname+"deltaPhi_det.png");

}

void Histos::Plot2D(bool debug) {

  // dump contents of the 2D histos
  if (debug) {
    Print2D(hsn_thmuVsthe_NLO_ref);
    Print2D(hsn_thL_vs_thR_NLO_ref);
    Print2D(hsn_thmuVsthe_NLO_lep);
    Print2D(hsn_thmuVsthe_NLO_templ[0]);
    Print2D(hsn_thmuVsthe_NLO_templ[440]);
    Print2D(hsn_thmuVsthe_LO_lep);
    Print2D(hsn_thmuVsthe_LO_templ[0]);
    Print2D(hsn_thmuVsthe_LO_templ[440]);
  }
  
  gStyle->SetOptStat(11);
  TString pname = dirname+"/"+prefix;

  TCanvas cz_thmuVsthe("cz_thmuVsthe", "Muon vs Electron angle", 700, 500);
  hsn_thmuVsthe_NLO_ref->Draw("COLZ");
  cz_thmuVsthe.SaveAs(pname+"thmuVsthe_det_colz.png");

  TCanvas cz_thL_vs_thR("cz_thL_vs_thR", "Theta_L vs Theta_R angle", 700, 700);
  hsn_thL_vs_thR_NLO_ref->Draw("COLZ");
  cz_thL_vs_thR.SaveAs(pname+"thL_vs_thR_det_colz.png");
}

void Histos::Print2D(TH2D *h) const {

  cout<<endl<<"Dumping contents of histogram: "<<h->GetName()<<endl;
  Int_t nx = h->GetNbinsX();
  Int_t ny = h->GetNbinsY();
  cout<<" nr of X bins = "<<nx <<", ny of Y bins = "<<ny <<endl;
  cout<<"  ix " <<", iy  :"<< " content, error " <<endl;
  for (Int_t ix=1; ix<nx+1; ix++) {
    for (Int_t iy=1; iy<ny+1; iy++) {
      double c = h->GetBinContent(ix,iy);
      double e = h->GetBinError(ix,iy);
      if (c<0) cout<<"  "<<ix <<" , "<<iy <<"  : "<< c << ",  "<< e << "*** NEGATIVE ***" <<endl; 
      else     cout<<"  "<<ix <<" , "<<iy <<"  : "<< c << ",  "<< e <<endl;
    }
  }
}

void Histos::Ratio_Resolution() {
  hsratio_t13_Resolution->Divide(hsn_t13, hn_t13);
  hsratio_t24_Resolution->Divide(hsn_t24, hn_t24);
  hsratio_the_Resolution->Divide(hsn_the, hn_the);
  hsratio_thmu_Resolution->Divide(hsn_thmu, hn_thmu);
  hsratio_x13_Resolution->Divide(hsn_x13, hn_x13);
  hsratio_x24_Resolution->Divide(hsn_x24, hn_x24);
  
  hsratio_t13_Resolution_LO->Divide(hsn_t13_LO, hn_t13_LO);
  hsratio_t24_Resolution_LO->Divide(hsn_t24_LO, hn_t24_LO);
  hsratio_the_Resolution_LO->Divide(hsn_the_LO, hn_the_LO);
  hsratio_thmu_Resolution_LO->Divide(hsn_thmu_LO, hn_thmu_LO);
  hsratio_x13_Resolution_LO->Divide(hsn_x13_LO, hn_x13_LO);
  hsratio_x24_Resolution_LO->Divide(hsn_x24_LO, hn_x24_LO);
}

void Histos::Ratio_FullRunning() {
  hratio_t13_FullRunning->Divide(hn_t13, hn_t13_norun);
  hratio_t24_FullRunning->Divide(hn_t24, hn_t24_norun);
  hratio_Ee_FullRunning->Divide(hn_Ee, hn_Ee_norun);
  hratio_Emu_FullRunning->Divide(hn_Emu, hn_Emu_norun);
  hratio_the_FullRunning->Divide(hn_the, hn_the_norun);
  hratio_thmu_FullRunning->Divide(hn_thmu, hn_thmu_norun);
  hratio_x13_FullRunning->Divide(hn_x13, hn_x13_norun);
  hratio_x24_FullRunning->Divide(hn_x24, hn_x24_norun);
  
  hsratio_t13_FullRunning->Divide(hsn_t13, hsn_t13_norun);
  hsratio_t24_FullRunning->Divide(hsn_t24, hsn_t24_norun);
  hsratio_Ee_FullRunning->Divide(hsn_Ee, hsn_Ee_norun);
  hsratio_Emu_FullRunning->Divide(hsn_Emu, hsn_Emu_norun);
  hsratio_the_FullRunning->Divide(hsn_the, hsn_the_norun);
  hsratio_thmu_FullRunning->Divide(hsn_thmu, hsn_thmu_norun);
  hsratio_x13_FullRunning->Divide(hsn_x13, hsn_x13_norun);
  hsratio_x24_FullRunning->Divide(hsn_x24, hsn_x24_norun);
}

void Histos::Ratio_HadronicRunning() {
  hratio_t13_HadronicRunning->Divide(hn_t13, hn_t13_lep);
  hratio_t24_HadronicRunning->Divide(hn_t24, hn_t24_lep);
  hratio_Ee_HadronicRunning->Divide(hn_Ee, hn_Ee_lep);
  hratio_Emu_HadronicRunning->Divide(hn_Emu, hn_Emu_lep);
  hratio_the_HadronicRunning->Divide(hn_the, hn_the_lep);
  hratio_thmu_HadronicRunning->Divide(hn_thmu, hn_thmu_lep);
  hratio_x13_HadronicRunning->Divide(hn_x13, hn_x13_lep);
  hratio_x24_HadronicRunning->Divide(hn_x24, hn_x24_lep);
  
  hsratio_t13_HadronicRunning->Divide(hsn_t13, hsn_t13_lep);
  hsratio_t24_HadronicRunning->Divide(hsn_t24, hsn_t24_lep);
  hsratio_Ee_HadronicRunning->Divide(hsn_Ee, hsn_Ee_lep);
  hsratio_Emu_HadronicRunning->Divide(hsn_Emu, hsn_Emu_lep);
  hsratio_the_HadronicRunning->Divide(hsn_the, hsn_the_lep);
  hsratio_thmu_HadronicRunning->Divide(hsn_thmu, hsn_thmu_lep);
  hsratio_x13_HadronicRunning->Divide(hsn_x13, hsn_x13_lep);
  hsratio_x24_HadronicRunning->Divide(hsn_x24, hsn_x24_lep);
}

void Histos::Do_Ratios() {
  Ratio_Resolution();
  Ratio_FullRunning();
  Ratio_HadronicRunning();    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following functions have been used for specific tasks
//
void Histos::RatioFinal(Long_t n_events, TFile *fp, TString hname) {
  
  cout<<"\n*** Computing Final ratios (with Projected statistics) for histogram: "<< hname <<endl;

  Long_t nevtot = n_events > 0 ? n_events : mcsums.Nevgen;
  Double_t wnorm = UNWGT ? mcsums.Xsec : Wnorm;
  Double_t norm = wnorm / static_cast<Double_t>(nevtot);

  cout<< "\n nevtot = mcsums.Nevgen = "<<nevtot<<endl;
  cout<< " wnorm = "<<wnorm<<endl;
  cout<< " norm = wnorm / nevtot = "<< norm << endl;

  cout<< " Reference Xsec = "<< mcsums.Xsec << " +/- " << mcsums.XsecErr << endl; 

  Double_t IntLumi = 1.5e10; // in microbarns^-1 // = 1.5*10^7 nb^-1
  cout<< " Projected Integrated Luminosity (paper) = "<< IntLumi <<" ub^-1"<<endl;

  Double_t Nnorm = norm*IntLumi;
  cout<< "final norm. factor to Countings for final Lumi = "<<Nnorm <<endl;

  // Ee>0.2 GeV
  //  Double_t SigmaRef = 1337.26;     // in microbarns
  //  Double_t SigmaRef_Err = 0.0172;
  // Ee>5 GeV
  //  Double_t SigmaRef = 45.9925;     // in microbarns
  //  Double_t SigmaRef_Err = 0.000486;
  //   cout<<"Reference Sigma : "<<SigmaRef<<" +- "<<SigmaRef_Err<<" microbarn"<<endl;

  cerr << "\n" << "Loading Reference root histograms from the existing file: "<< fp->GetName() <<endl;

  TH1D *hn = (TH1D*)fp->Get(hname)->Clone(hname+"_ref");
  Int_t nbinsm = hn->GetNbinsX();
  cout<<"n bins           : "<<nbinsm<<endl;
  cout<<"entries          : "<<hn->GetEntries()<<endl;
  cout<<"effective entries: "<<hn->GetEffectiveEntries()<<endl;
  cout<<"Reference integrated luminosity: "<< hn->GetEffectiveEntries() / mcsums.Xsec <<endl;
  cout<<"sum of weights (senza und/ovf): "<<hn->GetSumOfWeights()<<endl;
  Double_t integral, err;
  integral = hn->IntegralAndError(1,100,err);
  cout<<"Bin 1-100: integral, error: "<<integral<<", "<<err<<endl; //<<", effective integral: "<<effIntegral<<endl;
  // normalizzo a norm=wnorm/nevtot
  Double_t n_integral = integral*norm;
  Double_t n_err = err*norm;
  cout<<"(normalized) integr., err : "<<n_integral<<", "<<n_err<<endl;
  // ... add und/ovf
  cout<<"underflow content and error: "<<hn->GetBinContent(0)<<" +/- "<<hn->GetBinError(0)<<endl;
  cout<<"overflow  content and error: "<<hn->GetBinContent(nbinsm+1)<<" +/- "<<hn->GetBinError(nbinsm+1)<<endl;
  integral = hn->IntegralAndError(0,nbinsm+1,err);
  cout<<"Bin 0-101: integral, error: "<<integral<<", "<<err<<endl; //", effective integral: "<<effIntegral<<endl;
  n_integral = integral*norm;
  n_err = err*norm;
  cout<<"(normalized) integr., err : "<<n_integral<<", "<<n_err<<endl;

  // build the errors
  //  const TArrayD* sw2 = href_nthmu->GetSumw2();  // sono uguali al quadrato di GetBinError(im)

  TH1D *hrefRelErr = (TH1D*)hn->Clone(hname+"_RelErr");
  hrefRelErr->Reset();
  for (Int_t im=0; im<=nbinsm+1; ++im) {
    Double_t l_im = hn->GetBinLowEdge(im);
    Double_t n_im = hn->GetBinContent(im);
    Double_t err_im = hn->GetBinError(im);
    //Double_t e2_im = sw2->At(im);
    //Double_t neff = n_im*n_im /err_im*err_im;
    //Double_t reff = n_im/e2_im; // = neff/n_im 
    Double_t relerr_im = n_im > 0 ? err_im / n_im : 0; 
    Double_t neff_im = n_im > 0 ? 1./(relerr_im*relerr_im) : 0;
    Double_t reff_im = n_im > 0 ? neff_im / n_im : 0;

    hrefRelErr->SetBinContent(im, relerr_im);

    cout<<"bin "<<im<<" ("<<l_im<<")"<<", RelErr: "<<relerr_im<<"\n\t content: "<<n_im<<", error: "<<err_im 
      	<<", N effective entries: "<<neff_im<<", Ratio Eff/Gen: "<<reff_im<<endl;
  }

  // create histos with final number of countings
  TH1D *hn_final = (TH1D*)hn->Clone(hname+"_final");
  hn_final->Scale(Nnorm);

  TH1D *hfinalRelErr = (TH1D*)hn->Clone(hname+"_RelErr_final");
  hfinalRelErr->Reset();
  // build the projected errors
  for (Int_t im=0; im<=nbinsm+1; ++im) {
    Double_t n_im = hn_final->GetBinContent(im);
    Double_t relerr_im = n_im > 0 ? 1./sqrt(n_im) : 0; 
    hfinalRelErr->SetBinContent(im, relerr_im);
  }

  // rapporti tra gli errori^2 degli eventi generati e della proiezione
  cout<<"\n Ratios of errors for gen events and the projected full stat:"<<endl;
  for (Int_t im=0; im<=nbinsm+1; ++im) {
    Double_t relerr_gen  = hrefRelErr->GetBinContent(im);
    Double_t relerr_proj = hfinalRelErr->GetBinContent(im);
    Double_t ratio = relerr_gen / relerr_proj;
    Double_t fstat = ratio*ratio;
    cout<<"bin "<<im<<", GEN error: "<<relerr_gen<<", PROJ error: "<<relerr_proj<<", stat.factor = "<<fstat<<endl;
  }

  TString hname_ratioFull = hname;
  hname_ratioFull = (hname_ratioFull.Replace(0,2,"hratio"))+"_FullRunning";
  TString hname_ratioHad = hname;
  hname_ratioHad = (hname_ratioHad.Replace(0,2,"hratio"))+"_HadronicRunning";

  TH1D* h_ratio_FullRunning = (TH1D*)fp->Get(hname_ratioFull)->Clone(hname_ratioFull+"_ref");
  TH1D* h_ratio_HadronicRunning = (TH1D*)fp->Get(hname_ratioHad)->Clone(hname_ratioHad+"_ref");

  TH1D* hratioFinal_FullRunning = (TH1D*)h_ratio_FullRunning->Clone(hname_ratioFull+"_final");
  TH1D* hratioFinal_HadronicRunning =(TH1D*)h_ratio_HadronicRunning->Clone(hname_ratioHad+"_final");

  cout <<"\n test hratioFinal_HadronicRunning (before)"<<endl;
  if (hratioFinal_HadronicRunning->GetNbinsX() != nbinsm) {cout<<"***ERROR*** 2"<<endl; exit(2);}
  for (Int_t im=0; im<=nbinsm+1; ++im) {
    Double_t x_im   = hratioFinal_HadronicRunning->GetBinLowEdge(im);
    Double_t val_im = hratioFinal_HadronicRunning->GetBinContent(im);
    cout<<"bin "<<im<<" ("<<x_im<<"), value: "<<val_im<<endl;
  }
  // now reset the errors
  for (Int_t im=0; im<=nbinsm+1; ++im) {
    Double_t valf = hratioFinal_FullRunning->GetBinContent(im);
    Double_t errf = valf * hfinalRelErr->GetBinContent(im);
    valf += gRandom->Gaus(0.,errf);
    hratioFinal_FullRunning->SetBinContent(im, valf);
    hratioFinal_FullRunning->SetBinError(im, errf);

    Double_t valh = hratioFinal_HadronicRunning->GetBinContent(im);
    Double_t errh = valh * hfinalRelErr->GetBinContent(im);
    valh += gRandom->Gaus(0.,errh);
    hratioFinal_HadronicRunning->SetBinContent(im, valh);
    hratioFinal_HadronicRunning->SetBinError(im, errh);
  }
  cout <<"\n test hratioFinal_HadronicRunning (after error rescaling)"<<endl;
  if (hratioFinal_HadronicRunning->GetNbinsX() != nbinsm) {cout<<"***ERROR*** 3"<<endl; exit(3);}
  for (Int_t im=0; im<=nbinsm+1; ++im) {
    Double_t x_im   = hratioFinal_HadronicRunning->GetBinLowEdge(im);
    Double_t val_im = hratioFinal_HadronicRunning->GetBinContent(im);
    Double_t err_im = hratioFinal_HadronicRunning->GetBinError(im);
    Double_t relerr_im = val_im > 0 ? err_im / val_im : 0;
    cout<<"bin "<<im<<" ("<<x_im<<")"<<", value: "<<val_im<<", error: "<<err_im<<", RelErr: "<<relerr_im<<endl;
  }
}


void Histos::LoadCarlos(const vector<pair<string,Long_t> > & mc_inputs) {
  string line;
  istringstream stream;
  Double_t c_lowedge, c_content, c_error, c_dump;
  
  for (const auto & kdist : kine_distributions) {
    string name, carlosName;
    std::tie(name, carlosName, std::ignore) = kdist;
    string fname = carlosName+"_oal_100.txt";
    cerr << "Loading Carlos distribution: "<< fname << endl;
    
    Double_t xbinlow[100]{}, y[100]{}, ey[100]{};
    Long_t N_tot_weights = 0;
    
    for (const auto & idir : mc_inputs) {
      ifstream carlos_file(idir.first+fname);
      N_tot_weights += idir.second;
      Double_t fstat = static_cast<Double_t>(idir.second);
      Int_t ibin = 0;
      while (getline(carlos_file, line)) {
	if (ibin>=100) {
	  cerr << endl << "Error: more than 100 lines in input Carlos file: " << fname << endl;
	  exit(300);
	}
	stream = istringstream(line);
	//      stream >> c_lowedge >> c_content >> c_error >> c_dump;   // error computed a-la-MC (Carlo)
	stream >> c_lowedge >> c_content >> c_dump >> c_error;  // error computed as in root
	xbinlow[ibin] = c_lowedge;
	if (c_error > 0) {
	  y[ibin]       = y[ibin] + c_content*fstat;
	  ey[ibin]      = ey[ibin] + (c_error*fstat)*(c_error*fstat);
	}
	ibin++;
      }
      if (ibin<100) {
	cerr << endl << "Error: less than 100 lines in input Carlos file: " << fname << endl;
	exit(400);
      }
    }
    // combine all the samples and average
    Double_t f_ntotw = static_cast<Double_t>(N_tot_weights);
    for (Int_t jbin=0; jbin<100; ++jbin) {
      if (ey[jbin] >0) {
	y[jbin]=y[jbin]/f_ntotw;
	ey[jbin]=sqrt(ey[jbin])/f_ntotw;
      } else {y[jbin]=0; ey[jbin]=0;}
    }
    
    ofstream sumdist;
    string fnameSum = carlosName+"_total.txt";
    sumdist.open(fnameSum);
    for (Int_t kbin=0; kbin<100; ++kbin) {
      sumdist << setw(25) << setprecision(17) << xbinlow[kbin] 
	      << setw(25) << setprecision(17) << y[kbin]
	      << setw(25) << setprecision(17) << ey[kbin]
	      << endl;
    }
    sumdist.close();
  }
}

void Histos::LoadExtHistos(TFile * fp) {
  cerr << "\n" << "Loading root histograms from the existing file: "<<fp->GetName()<<endl;
  fp->GetObject("h_t13", h_t13);
  fp->GetObject("h_t24", h_t24);
  fp->GetObject("h_Ee", h_Ee);
  fp->GetObject("h_the", h_the);
  fp->GetObject("h_Emu", h_Emu);
  fp->GetObject("h_thmu", h_thmu);
  fp->GetObject("h_egamma", h_egamma);
  fp->GetObject("h_thgamma", h_thgamma);
}

void Histos::CompareWithCarlos(bool debug) {
  cerr << "\n" << "Comparing Carlo's distributions with my root histograms." << endl;
  
  for (const auto & kdist : kine_distributions) {
    string name, carlosName;
    TH1D** addr_my_hist;
    std::tie(name, carlosName, addr_my_hist) = kdist;
    TH1D* my_hist = *addr_my_hist;
    string fname = carlosName+"_total.txt";
    cout << "\n" << "comparing Carlos distribution : " << fname
	      << " with my histogram "<< my_hist->GetName() << endl;
    
    string line;
    istringstream stream;
    ifstream carlos_file(fname);
    Int_t i = 0;
    Double_t c_lowedge, c_content, c_error;
    Double_t lowedge, content, error;
    
    Double_t maxdif = 0.;
    Double_t maxdif_err = 0.;
    Int_t ibmax = 0;
    Int_t ibmaxe = 0;
    
    TGraph g_gio(100), g_gio_err(100), g_carlo(100), g_carlo_err(100), g_diff(100), g_diff_err(100);
    
    while (getline(carlos_file, line)) {
      i++;
      if (i>100) {
	cerr << endl << "Error: more than 100 lines in input Carlos file: " << fname << endl;
	return;
      }
      stream = istringstream(line);
      stream >> c_lowedge >> c_content >> c_error;
      
      g_carlo.SetPoint(i, c_lowedge, c_content);
      g_carlo_err.SetPoint(i, c_lowedge, c_error);
      
      lowedge = my_hist->GetBinLowEdge(i);
      content = my_hist->GetBinContent(i);
      error = my_hist->GetBinError(i);
      
      g_gio.SetPoint(i, lowedge, content);
      g_gio_err.SetPoint(i, lowedge, error);
      
      Double_t diff = c_content != 0. ? (content - c_content) /c_content : +1. ;
      if (c_content == 0. && content == 0.) diff = 0.;
      if (std::abs(diff) > maxdif) {
	ibmax = i;
	maxdif = std::abs(diff); 
      }
      
      Double_t diff_err = c_error != 0. ? (error - c_error) /c_error : +1. ;
      if (c_error == 0. && error == 0.) diff_err = 0.;
      if (std::abs(diff_err) > maxdif_err) {
	ibmaxe = i;
	maxdif_err = std::abs(diff_err); 
      }
      
      if (debug)
	cout<< "c_content, content = " << c_content << " " << content << " rel diff cont = "<< diff
		 << "\t c_error, error     = " << c_error << " " << error << " rel diff error = "<< diff_err  <<endl;
      
      g_diff.SetPoint(i, c_lowedge, diff);
      g_diff_err.SetPoint(i, c_lowedge, diff_err);
    }
    
    if (debug)
      cout << "*** MAX Relative Difference in contents : "<<  maxdif << ", at bin : "<< ibmax
		<< "*** MAX Relative Difference in errors : "<<  maxdif_err << ", at bin : "<< ibmaxe <<endl;
    
    TCanvas cc (carlosName.c_str(), carlosName.c_str());
    g_diff.SetTitle(("relative difference Carlo-root on: "+carlosName).c_str());
    g_diff.Draw("ALP");
    cc.SaveAs(("reldif_"+carlosName+".png").c_str());
    
    TCanvas cce ((carlosName+"_err").c_str(), (carlosName+"_err").c_str());
    g_diff_err.Draw("ALP");
    g_diff_err.SetTitle(("relative difference Carlo-root on: " +carlosName+" (error)").c_str());
    cce.SaveAs(("reldif_"+carlosName+"_err.png").c_str());
  }
}
