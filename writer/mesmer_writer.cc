/////////////////////////////////////////////////////////////////////////////////////////
//
// Generate MESMER MuE MC events and write them out in root format
//
//   Input parameters: filenames for
//   1: input MESMER card file 
//   2: root output
//
// G.Abbiendi  6/Oct/2021 previous writer code
// G.A.+C.Carloni Calame 4/Mar/2022: embedded MESMER generation 
/////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TTree.h"
#include "MuEtree.h"
#include "Mesmer.h"

using namespace std;

// set it to true to debug
bool DEBUG = false;
//bool DEBUG = true;

int main(int argc, char* argv[]) {

  if (argc != 3) {
    cout << "Usage : "
	 << argv[0] << " input_mesmer_card output_ROOT_filename \n";
    exit(100);
  }

  TString ofile(argv[2]);
  TFile *output_file = new TFile(ofile,"RECREATE");
  
  cout<<endl<<"Reading Mesmer run parameters from input file: "<< argv[1] <<endl;
  cout      <<"Writing events to output root file           : "<< argv[2] <<endl;

  MuE::Setup setup;

  // Initialise MESMER
  setup.MCpargen.InitandSetRunParams_mesmer(argv[1]);
  
  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("Root Application", &argc, argv);

  Int_t splitBranches = 2;
  
  TTree *partree = new TTree("MuEsetup","MuE MC parameters");
  partree->Branch("MuEparams",&setup,64000,splitBranches);

  MuE::Event event;
  TTree *tree = new TTree("MuEtree","MuE MC tree");
  tree->Branch("MuE",&event,64000,splitBranches);
  tree->SetAutoSave(500000);

  // number of requested events
  const Long_t HowManyEvts = setup.MCpargen.Nevreq;

  // Event generation loop 
  cout << "Start generating MESMER events" <<endl;
  Long_t nev = 0;
  while (nev < HowManyEvts) {

    double pmu[4];
    IncomingMuonMomentum_mesmer(pmu); // gets the initial state muon 4-momentum (E,px,py,pz) components
                                      // can be replaced by any 4-momentum "provider" function
    
    int ierr = event.GenerateEvent_mesmer(pmu);

    tree->Fill();

    if (ierr == 0) nev++;
    //    else event.Print();

    if (DEBUG) event.Print();
  }

  cout << "End generating MESMER events." <<endl;
  if (DEBUG) tree->Print();

  // Fill summary MC statistics
  setup.MCsums.SetEndofRun_mesmer();
  
  // number of generated weights (including zero-weight events)
  Long64_t Nentries = tree->GetEntries();
  if (Nentries != setup.MCsums.Nwgt) {
    cout << "\n *** ERROR *** mismatch: number of entries in the tree is "<< Nentries 
	 << ", but the number of generated weights is "<<setup.MCsums.Nwgt << endl;
  }

  setup.MCsums.Print();
  
  // fill the parameter tree (1 entry per mc run)
  partree->Fill();
  if (DEBUG) partree->Print();

  output_file->Write();
  output_file->Close();
  cout << "Generated " << HowManyEvts << " events, written to file: "<< ofile << endl;

  if (!gROOT->IsBatch()) app->Run();

  return 0;
}
