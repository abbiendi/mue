#ifndef Histos_H
#define Histos_H

////////////////////////////////////////////////////////////////
// class for MuE histograms
//
// G.Abbiendi  5/Jul/2018 
////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TArrayD.h"
#include "MuEtree.h" 
#include "MuEana.h"
#include "Inputs.h"

namespace MuE {
  
  class Histos {
    
  public:
    Histos(const MCpara & pargen, const AN_Input & paran, const TString & nameDir="MuE", const TString & namePrefix="", const bool make_dirs=true, bool debug=false);
    virtual ~Histos(){};
    void Fill(const Event & event, const MuEana & a);
    void FillTemplates(const Event & event, const MuEana & a);
    void FillTemplatesLO(const MuEana & a);
    inline void SetSums(const MCstat & s) {mcsums = s;}
    void Normalize();
    void FillResolution(const MuEana & a, bool debug=false);

    void Fit();
    void Plot();
    void Plot2D(bool debug=false);
    void PlotResolutions();
    void Ratio_Resolution();
    void Ratio_FullRunning();
    void Ratio_HadronicRunning();
    void Do_Ratios();
    void RatioFinal(Long_t n_events, TFile *fp, TString hname);
    void LoadExtHistos(TFile *fp);
    void LoadCarlos(const std::vector<std::pair<std::string,Long_t> > & mc_inputs);
    void CompareWithCarlos(bool debug=false);

  private:
    void hdraw(TH1D* hist);
    void Print2D(TH2D *h) const;

    typedef std::tuple<std::string,std::string,TH1D**> kine_dist; //name, carlosName, my_hist
    std::vector<kine_dist> kine_distributions;

    // parameters of the MC Generator (extracted from pargen)
    const Double_t & Emin_e;    // minimum electron energy in Carlo's distributions
    const Double_t & charge_mu; // muon charge
    const Double_t & mass_mu;   // muon mass
    const Double_t & mass_e;    // electron mass
    const Double_t & invalfa0;  // 1/alpha(0)
    const Bool_t & UNWGT;       // true for unweighted events
    const Double_t & Wnorm;     // normalizing factor in xsection 
    MCstat mcsums;              // Sums of the MC Generator
    //
    // parameters for analysis (extracted from paran)
    const bool & angleCorrPlots_finebins; // theta_mu vs theta_e in fine bins
    const bool & doTemplates;        // produce 2D template histograms
    const double & thetaMax;         // max theta for event selection (geometric acceptance)
    const double & theMaxTemp;       // max theta for template histograms
    const int & iLumi;               // select Int.Lumi (0:Nominal; 1:LowLumi(TR2021))
    const unsigned int & rangeSigma; // range around exp. ave (as number of +/-sigmas on each axis)
    const unsigned int & divSigma;   // number of divisions within one sigma interval
    const int & iparam;              // Da_had parameterization:  0:pol2 ; 1:LL ; 2:LLmod
    const bool & doEnergyScale;      // produce plots for E scale calibration
    //
    //--- calculated from input parameters
    unsigned int ngrid;
    double Kref, Mref, KerrRef, MerrRef; // reference expected central values and errors
    double dKu, dMu;
    //---

    const TString dirname;          // directory name where MuE results are going
    const TString prefix;            // prefix in histogram/canvas/plot names
    const bool makedirs;

    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////// CROSS SECTIONS (with Ideal Variables and Full Running) /////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    TH1D* h_t13;
    TH1D* h_t24;
    TH1D* h_Ee;
    TH1D* h_Emu;
    TH1D* h_the;
    TH1D* h_thmu;
    TH1D* h_egamma;
    TH1D* h_thgamma;
    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// EVENT DISTRIBUTIONS (weighted countings) /////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    //
    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// FULL RUNNING (Ideal Variables) /////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    // number of events (sum of weights with Full Running of alpha)
    TH1D* hn_t13;
    TH1D* hn_t24;
    TH1D* hn_Ee;
    TH1D* hn_Emu;
    TH1D* hn_the;
    TH1D* hn_thmu;
    TH1D* hn_egamma;
    TH1D* hn_thgamma;
    TH1D* hn_egamma_CoM;
    TH1D* hn_x13;
    TH1D* hn_x24;
    TH1D* hn_deltaPhi;
    TH2D* hn_thmu_vs_the;
    TH1D* hn_openingAngle;
    TH1D* hn_tripleProduct;
    //////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// FULL RUNNING (Smeared Variables) /////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // Histograms with Detector smearing
    // number of events (sum of weights with Full Running of alpha)
    TH1D* hsn_t13;
    TH1D* hsn_t24;
    TH1D* hsn_Ee;
    TH1D* hsn_Emu;
    TH1D* hsn_the;
    TH1D* hsn_thmu;
    TH1D* hsn_egamma;
    TH1D* hsn_thgamma;
    TH1D* hsn_egamma_CoM;
    TH1D* hsn_x13;
    TH1D* hsn_x24;
    TH1D* hsn_deltaPhi;
    TH2D* hsn_thmu_vs_the;
    TH1D* hsn_openingAngle;
    TH1D* hsn_tripleProduct;
    //
    //TH1D* hsn_theX;
    //////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// NO RUNNING (Ideal variables) /////////////
    //////////////////////////////////////////////////////////////////////////////////
    // number of events (sum of weights with NO Running of alpha)
    TH1D* hn_t13_norun;
    TH1D* hn_t24_norun;
    TH1D* hn_Ee_norun;
    TH1D* hn_Emu_norun;
    TH1D* hn_the_norun;
    TH1D* hn_thmu_norun;
    TH1D* hn_x13_norun;
    TH1D* hn_x24_norun;
    ///////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// NO RUNNING (Smeared variables) ////////////
    //////////////////////////////////////////////////////////////////////////////////
    // number of events (sum of weights with NO Running of alpha)
    TH1D* hsn_t13_norun;
    TH1D* hsn_t24_norun;
    TH1D* hsn_Ee_norun;
    TH1D* hsn_Emu_norun;
    TH1D* hsn_the_norun;
    TH1D* hsn_thmu_norun;
    TH1D* hsn_x13_norun;
    TH1D* hsn_x24_norun;

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// LEPTONIC RUNNING (Ideal variables) //////////
    /////////////////////////////////////////////////////////////////////////////////////
    // number of events (sum of weights with only Leptonic Running of alpha)
    TH1D* hn_t13_lep;
    TH1D* hn_t24_lep;
    TH1D* hn_Ee_lep;
    TH1D* hn_Emu_lep;
    TH1D* hn_the_lep;
    TH1D* hn_thmu_lep;
    TH1D* hn_x13_lep;
    TH1D* hn_x24_lep;
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// LEPTONIC RUNNING (Smeared variables) ////////
    /////////////////////////////////////////////////////////////////////////////////////
    // number of events (sum of weights with only Leptonic Running of alpha)
    TH1D* hsn_t13_lep;
    TH1D* hsn_t24_lep;
    TH1D* hsn_Ee_lep;
    TH1D* hsn_Emu_lep;
    TH1D* hsn_the_lep;
    TH1D* hsn_thmu_lep;
    TH1D* hsn_x13_lep;
    TH1D* hsn_x24_lep;

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// LEADING ORDER (Ideal variables) //////////
    /////////////////////////////////////////////////////////////////////////////////////
    // number of events (sum of weights for LO)
    TH1D* hn_t13_LO;
    TH1D* hn_t24_LO;
    TH1D* hn_Ee_LO;
    TH1D* hn_Emu_LO;
    TH1D* hn_the_LO;
    TH1D* hn_thmu_LO;
    TH2D* hn_thmu_vs_the_LO;
    TH1D* hn_x13_LO;
    TH1D* hn_x24_LO;
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// LEADING ORDER (Smeared variables) ///////////
    /////////////////////////////////////////////////////////////////////////////////////
    // number of events (sum of weights for LO)
    TH1D* hsn_t13_LO;
    TH1D* hsn_t24_LO;
    TH1D* hsn_Ee_LO;
    TH1D* hsn_Emu_LO;
    TH1D* hsn_the_LO;
    TH1D* hsn_thmu_LO;
    TH2D* hsn_thmu_vs_the_LO;
    TH1D* hsn_x13_LO;
    TH1D* hsn_x24_LO;
    //
    //    TH1D* hsn_theX_LO;
    ////////////////////////////////////////////////////////////////////////////////////
    // Angular Resolution plots
    ////////////////////////////////////////////////////////////////////////////////////
    TH1D* h_resol_thex;
    TH1D* h_resol_they;
    TH1D* h_resol_theSpaceXY;
    TH1D* h_resol_theSpace;
    TH1D* h_resol_thex5;
    TH1D* h_resol_they5;
    TH1D* h_resol_theSpace5XY;
    TH1D* h_resol_theSpace5;
    TH1D* h_resol_thmux;
    TH1D* h_resol_thmuy;
    TH1D* h_resol_thmuSpaceXY;
    //
    TH2D* h_resol_thexVsE;
    TH2D* h_resol_theyVsE;
    TH2D* h_resol_thmuxVsE;
    TH2D* h_resol_thmuyVsE;
    //
    TH1D* h_resol_thexVsE_sigma;
    TH1D* h_resol_theyVsE_sigma;
    TH1D* h_resol_thmuxVsE_sigma;
    TH1D* h_resol_thmuyVsE_sigma;

    //================================================================================
    //=================== Ratios Smeared / Ideal Variables (with full running) =======
    //================================================================================
    TH1D* hsratio_t13_Resolution;
    TH1D* hsratio_t24_Resolution;
    TH1D* hsratio_the_Resolution;
    TH1D* hsratio_thmu_Resolution;
    TH1D* hsratio_x13_Resolution;
    TH1D* hsratio_x24_Resolution;
    //================================================================================
    //====================== Ratios Smeared / Ideal Variables (for LO) ===============
    //================================================================================
    TH1D* hsratio_t13_Resolution_LO;
    TH1D* hsratio_t24_Resolution_LO;
    TH1D* hsratio_the_Resolution_LO;
    TH1D* hsratio_thmu_Resolution_LO;
    TH1D* hsratio_x13_Resolution_LO;
    TH1D* hsratio_x24_Resolution_LO;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++ Ratios Run / NO Run = Full Running (Ideal variables) ++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TH1D* hratio_t13_FullRunning;
    TH1D* hratio_t24_FullRunning;
    TH1D* hratio_Ee_FullRunning;
    TH1D* hratio_Emu_FullRunning;
    TH1D* hratio_the_FullRunning;
    TH1D* hratio_thmu_FullRunning;
    TH1D* hratio_x13_FullRunning;
    TH1D* hratio_x24_FullRunning;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++ Ratios Run / NO Run = Full Running (Smeared variables) ++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TH1D* hsratio_t13_FullRunning;
    TH1D* hsratio_t24_FullRunning;
    TH1D* hsratio_Ee_FullRunning;
    TH1D* hsratio_Emu_FullRunning;
    TH1D* hsratio_the_FullRunning;
    TH1D* hsratio_thmu_FullRunning;
    TH1D* hsratio_x13_FullRunning;
    TH1D* hsratio_x24_FullRunning;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++ Ratios Run / Leptonic Run = Hadronic Running (Ideal variables) +++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TH1D* hratio_t13_HadronicRunning;
    TH1D* hratio_t24_HadronicRunning;
    TH1D* hratio_Ee_HadronicRunning;
    TH1D* hratio_Emu_HadronicRunning;
    TH1D* hratio_the_HadronicRunning;
    TH1D* hratio_thmu_HadronicRunning;
    TH1D* hratio_x13_HadronicRunning;
    TH1D* hratio_x24_HadronicRunning;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++ Ratios Run / Leptonic Run = Hadronic Running (Smeared variables) ++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TH1D* hsratio_t13_HadronicRunning;
    TH1D* hsratio_t24_HadronicRunning;
    TH1D* hsratio_Ee_HadronicRunning;
    TH1D* hsratio_Emu_HadronicRunning;
    TH1D* hsratio_the_HadronicRunning;
    TH1D* hsratio_thmu_HadronicRunning;
    TH1D* hsratio_x13_HadronicRunning;
    TH1D* hsratio_x24_HadronicRunning;

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // TEMPLATE HISTOS
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    //
    //**************************************************************************
    // Next-To-Leading Order
    //**************************************************************************
    // NLO theta_muon vs theta_electron
    //
    // NLO distributions with full running (reference=Jegerlehner)
    TH2D *hn_thmuVsthe_NLO_ref;
    TH2D *hsn_thmuVsthe_NLO_ref;
    // NLO distributions with leptonic-only running
    TH2D *hn_thmuVsthe_NLO_lep;
    TH2D *hsn_thmuVsthe_NLO_lep;
    // TEMPLATES for parameterized hadronic running
    std::vector<TH2D* > hn_thmuVsthe_NLO_templ;
    std::vector<TH2D* > hsn_thmuVsthe_NLO_templ;
    //
    // NLO cross sections with full running (reference=Jegerlehner)
    TH2D *h_thmuVsthe_NLO_ref;
    TH2D *hs_thmuVsthe_NLO_ref;

    // NLO :  theta_Left vs theta_Right  (NO ID)
    //
    // NLO distributions with full running (reference=Jegerlehner)
    TH2D *hn_thL_vs_thR_NLO_ref;
    TH2D *hsn_thL_vs_thR_NLO_ref;
    // NLO distributions with leptonic-only running
    TH2D *hn_thL_vs_thR_NLO_lep;
    TH2D *hsn_thL_vs_thR_NLO_lep;
    // TEMPLATES for parameterized hadronic running
    std::vector<TH2D* > hn_thL_vs_thR_NLO_templ;
    std::vector<TH2D* > hsn_thL_vs_thR_NLO_templ;
    //
    // NLO cross sections with full running (reference=Jegerlehner)
    TH2D *h_thL_vs_thR_NLO_ref;
    TH2D *hs_thL_vs_thR_NLO_ref;

    //**************************************************************************
    // Leading Order
    //**************************************************************************
    // LO theta_muon vs theta_electron
    //
    // LO distributions with full running (reference=Jegerlehner)
    TH2D *hn_thmuVsthe_LO_ref;
    TH2D *hsn_thmuVsthe_LO_ref;
    // LO distributions with leptonic-only running
    TH2D *hn_thmuVsthe_LO_lep;
    TH2D *hsn_thmuVsthe_LO_lep;    
    // TEMPLATES for parameterized hadronic running
    std::vector<TH2D* > hn_thmuVsthe_LO_templ;
    std::vector<TH2D* > hsn_thmuVsthe_LO_templ;

    // LO cross sections with full running (reference=Jegerlehner)
    TH2D *h_thmuVsthe_LO_ref;
    TH2D *hs_thmuVsthe_LO_ref;

    // LO :  theta_Left vs theta_Right  (NO ID)
    //
    // LO distributions with full running (reference=Jegerlehner)
    TH2D *hn_thL_vs_thR_LO_ref;
    TH2D *hsn_thL_vs_thR_LO_ref;
    // LO distributions with leptonic-only running
    TH2D *hn_thL_vs_thR_LO_lep;
    TH2D *hsn_thL_vs_thR_LO_lep;    
    // TEMPLATES for parameterized hadronic running
    std::vector<TH2D* > hn_thL_vs_thR_LO_templ;
    std::vector<TH2D* > hsn_thL_vs_thR_LO_templ;

    // LO cross sections with full running (reference=Jegerlehner)
    TH2D *h_thL_vs_thR_LO_ref;
    TH2D *hs_thL_vs_thR_LO_ref;

    ////////////////////////////////////////////////////////////////
    // Plots for Energy scale calibration (around equal-angle point)
    //
    TH1D *hsn_thA_ref;
    TH2D *hsn_thD_vs_thA_ref;
    TH2D *hsn_thDs_vs_thA_ref;
    TH1D *hsn_thA_ref_cut2;
    TH2D *hsn_thD_vs_thA_ref_cut2;
    TH2D *hsn_thDs_vs_thA_ref_cut2;
    TH1D *hsn_thA_LO_ref;
    TH2D *hsn_thD_vs_thA_LO_ref;
    TH2D *hsn_thDs_vs_thA_LO_ref;
    TH1D *hsn_thA_LO_ref_cut2;
    TH2D *hsn_thD_vs_thA_LO_ref_cut2;
    TH2D *hsn_thDs_vs_thA_LO_ref_cut2;
    //
    //=========================================================
    // Plots for Energy scale calibration (for chi2 scan)
    // NLO
    TH2D *hesn_thmu_vs_the;
    TH2D *hesn_thL_vs_thR;
    TH2D *hesn_thD_vs_thA;
    TH2D *hesn_thD_vs_alpha;
    TH2D *hesn_alpha_vs_thA;
    // LO
    TH2D *hesn_thmu_vs_the_LO;
    TH2D *hesn_thL_vs_thR_LO;
    TH2D *hesn_thD_vs_thA_LO;
    TH2D *hesn_thD_vs_alpha_LO;
    TH2D *hesn_alpha_vs_thA_LO;

  };
  
}

#endif
