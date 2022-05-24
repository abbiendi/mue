#include <string>
#include <iostream>
#include "Analysis.h"

using namespace MuE;
using namespace std;

void Analysis::BeginJob()
{
  string outhist = parmain.output_dir+"/results.root";
  output_hist_file = new TFile(outhist.data(),"RECREATE");
  
  if (paran.makeTree) 
    {
      string outtree = parmain.output_dir+"/outtree.root";
      output_tree_file = new TFile(outtree.data(),"RECREATE");
      atree = new TTree("atree","Analysis output tree");
      Int_t splitBranches = 2;
      atree->Branch("event",&myAna,64000,splitBranches);
    }
  
  output_hist_file->cd();
  
  cout<<"\n"<<"Analysis Inputs: selection: theta_e < "<< paran.thetaMax << " mrad"<<endl;
  cout<<"\n"<<"Analysis Inputs: MAX theta for template histos = "<< paran.theMaxTemp << " mrad"<<endl;

  histos = new Histos(pargen, paran, parmain.output_dir);

  cout<<"\n"<<"\t adding selection Ee > 20 GeV"<<endl;
  histos_e20 = new Histos(pargen, paran, parmain.output_dir, "e20_", false);

  cout<<"\n"<<"\t adding selection |Dphi|<40 mrad"<<endl;
  histos_dph40 = new Histos(pargen, paran, parmain.output_dir, "dph40_", false);
}


void Analysis::Analyze(const MuE::Event & event, const MuE::FastSim & fs)
{
  const PxPyPzEVector & p_mu_in = fs.Get_p_mu_in();
  const PxPyPzEVector & p_e_in = fs.Get_p_e_in();
  const PxPyPzEVector & p_mu_out = fs.Get_p_mu_out();
  const PxPyPzEVector & p_e_out = fs.Get_p_e_out();
  const std::vector<PxPyPzEVector> & V_photons = fs.Get_photons();

  const MuE::KineVars & genKin = fs.GetGenKin();
  const MuE::KineVars & detKin = fs.GetDetKin();
  const MuE::Photon & photon = fs.GetPhoton();

  // selection (on detector-level variables)
  if (!Select(detKin)) return;

  // filling my analysis variables
  myAna.RunNr = event.RunNr;
  myAna.EventNr = event.EventNr;
  myAna.wgt_full = event.wgt_full;
  myAna.wgt_norun = event.wgt_norun;
  myAna.wgt_lep = event.wgt_lep;
  myAna.wgt_LO = event.wgt_LO;
  myAna.wgt_NLO = event.wgt_NLO;
  myAna.p_mu_in = p_mu_in;     
  myAna.p_e_in = p_e_in;     
  myAna.p_mu_out = p_mu_out;     
  myAna.p_e_out = p_e_out;     
  myAna.V_photons = V_photons;     
  myAna.genKin = genKin;
  myAna.detKin = detKin;
  myAna.photon = photon;
  
  // filling my analysis tree
  if (paran.makeTree) atree->Fill();
  
  // filling my histos
  histos->Fill(event, myAna);
  //  
  // selections to study energy calibration
  if (detKin.Ee > 20.) histos_e20->Fill(event, myAna);
  if (std::abs(detKin.deltaPhi)<0.040) histos_dph40->Fill(event, myAna);
}


bool Analysis::Select(const KineVars & e)
{
  bool a = e.Ee > 1e-6 && e.the < paran.thetaMax;
  return a;
}


void Analysis::EndJob(const MCstat & mcsums, const FastSim & fs, std::vector<std::pair<std::string,Long_t> > *mc_inputs) 
{
  Long_t n_events = mcsums.Nevgen;

  if (n_events >0) {
    cerr << "Final histogram normalizations, ratios, fits, plots "<<endl<<endl;
    histos->SetSums(mcsums);
    histos->Normalize();
    histos->Do_Ratios();
    histos->Plot();
    if (paran.doTemplates) histos->Plot2D();
    //    if (paran.doTemplates) histos->Plot2D(true);
    histos->Fit();
    histos->PlotResolutions();

    // selections to study energy calibration
    histos_e20->SetSums(mcsums);
    histos_e20->Normalize();
    histos_e20->Do_Ratios();
    histos_e20->Plot();
    if (paran.doTemplates) histos_e20->Plot2D();
    histos_e20->Fit();
    histos_e20->PlotResolutions();

    histos_dph40->SetSums(mcsums);
    histos_dph40->Normalize();
    histos_dph40->Do_Ratios();
    histos_dph40->Plot();
    if (paran.doTemplates) histos_dph40->Plot2D();
    histos_dph40->Fit();
    histos_dph40->PlotResolutions();

    // write out root tree
    if (paran.makeTree) {
      atree->Print();
      output_tree_file->Write();
    }
    
    // write out hist root file
    output_hist_file->Write();
    
    cerr <<endl<< "Processed " << n_events << " events, histograms written to file: results.root" << endl;
  }

  else {
    // load histos from existing external file
    TFile *fp = new TFile(parmain.histo_ifname.c_str());

    //    output_hist_file->cd();
    
    TFile *projfile = new TFile("projstat.root","RECREATE");
    projfile->cd();
    histos->RatioFinal(n_events, fp, "hn_thmu");
    histos->RatioFinal(n_events, fp, "hn_the");
    projfile->Write();

    cerr <<"\n"<< "*** WARNING: missing code to normalize histograms to be used in analysis !!! \n"<<endl;

    // test numerical precision comparing results from FORTRAN code to ROOT
    if (mc_inputs != nullptr) {
      histos->LoadExtHistos(fp);
      histos->LoadCarlos(*mc_inputs);
      histos->CompareWithCarlos();
    }
  }

}
