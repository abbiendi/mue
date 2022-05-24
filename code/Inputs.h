#ifndef Inputs_H
#define Inputs_H

namespace MuE {

  class Inputs {
  public:
    Inputs(){};
    virtual ~Inputs(){};
    void configure(const std::string & cfi, bool debug=false);
    virtual void readcfi(std::ifstream & ifile, bool debug=false){};
    void readin(std::ifstream & ifile, auto & vvv, bool debug=false);
  };

  class One_Input : public Inputs {
  public:
    One_Input(){};
    virtual ~One_Input(){};
    void readcfi(std::ifstream & ifile, bool debug=false);

    std::string input_dirs_file; // file with list of dir. paths with input MC events
    long n_events;               // n.events to be processed (N>0:N; 0:all; <0: read histograms from existing results)
    std::string ifname;          // input string in event file name
    std::string histo_ifname;    // complete path of file with input histos (for n_events<0)
    std::string output_dir;      // output directory name 
    std::string fastSim_ifname;  // FastSim input cfg 
    std::string analysis_ifname; // Analysis input cfg
  };

  class FS_Input : public Inputs {
  public:
    FS_Input(){};    
    virtual ~FS_Input(){};
    void readcfi(std::ifstream & ifile, bool debug=false);

    int model;         // detector resolution model: 0=simplest-2par; 1=Antonio's 3 par
    int MSopt;         // 0:default; 1:MS only on X plane; 2: polar angle smearing
    bool twosteps;     // if true in model 0 apply only multiple scattering
    double thickness;  // total material thickness for model 0 (in X0) / 1.5cm Be = 0.0425 X0
    double resolution; // intr.ang.resol. for model 0 (in mrad) // optimal 10um/50cm=0.02 mrad
    double zBias;      // signed bias of station length (um)
    double zSigma;     // uncertainty of station length (um)
  };

  class AN_Input : public Inputs {
  public:
    AN_Input(){};   
    virtual ~AN_Input(){};
    void readcfi(std::ifstream & ifile, bool debug=false);
    inline unsigned int ngrid() const {return 2*rangeSigma*divSigma + 1;}
    void hadrpars(double & Kref, double & Mref, double & KerrRef, double & MerrRef, double & dKu, double & dMu) const;
  
    double thetaMax;         // max theta for event selection (geometric acceptance)
    double theMaxTemp;       // max theta for templates
    unsigned int rangeSigma; // range around exp. ave (as number of +/-sigmas on each axis)
    unsigned int divSigma;   // number of divisions within one sigma interval
    bool doTemplates;        // produce 2D template histograms
    int iLumi;               // select Int.Lumi (0:Nominal; 1:LowLumi(TR2021))
    int iparam;              // Da_had parameterization:  0:pol2 ; 1:LL ; 2:LLmod
    bool makeTree;           // make tree with event variables
    bool angleCorrPlots_finebins; // theta_mu vs theta_e in fine bins
    bool doEnergyScale;      // produce plots for E scale calibration
  };
  
}

#endif
