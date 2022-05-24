#ifndef MuEtree_H
#define MuEtree_H

///////////////////////////////////////////////
// Data Formats for MuE MC events
//
// G.Abbiendi  v1 5/Jul/2018
//             v2 6/Oct/2021
//             v3 18/Feb/2022
///////////////////////////////////////////////

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include "TROOT.h"

namespace MuE {
  
  class MCpara {

  public:
    Int_t SampleTag;
    std::string program_version;
    Int_t process_ID;
    std::string running_on;
    std::string start_time;
    Long_t Nevreq;
    Bool_t UNWGT;
    std::string Mode;
    Int_t nphotmode;
    Int_t radmuch, radelch;
    Int_t rnd_ext, rnd_int;
    Double_t Ebeam, EbeamRMS;
    Bool_t EXTBEAM;
    Double_t charge_mu;
    Double_t mass_mu, mass_e;
    Double_t invalfa0;
    Double_t Emin_e;
    Double_t thmin_e, thmax_e;
    Double_t thmin_mu, thmax_mu;
    Double_t Ethr, ththr;
    Int_t i_acopl;
    Double_t cut_acopl;
    Int_t i_elast;
    Double_t cut_elast;
    Double_t Wnorm; 
    Double_t Wmax;
    Bool_t READ_COEF;
    Int_t ihadVP;
    Int_t ihadVPfl;
    Int_t Nsearch;
    Int_t Ndistr;
    Int_t isync;
    Double_t eps; 
    Double_t phmass; 

    MCpara(std::ifstream &, bool debug=false);

    MCpara():
      SampleTag(0),process_ID(0),Nevreq(0),UNWGT(false),nphotmode(-1),radmuch(1),radelch(1),rnd_ext(0),rnd_int(0),
      Ebeam(0),EbeamRMS(0),EXTBEAM(false),charge_mu(0),mass_mu(0),mass_e(0),invalfa0(0),
      Emin_e(0),thmin_e(0),thmax_e(0),thmin_mu(0),thmax_mu(0),Ethr(0),ththr(0),i_acopl(0),cut_acopl(0),i_elast(0),cut_elast(0),
      Wnorm(0),Wmax(0),READ_COEF(true),ihadVP(0),ihadVPfl(0),Nsearch(0),Ndistr(0),isync(0),eps(0),phmass(0)
    {};

    void InitandSetRunParams_mesmer(char * input_mesmer);

    void Print() const;
    void Dump() const;

    virtual ~MCpara(){};
    
    ClassDef(MCpara,3)
  };

  class MCstat {

  public:
    Long_t Nevgen;
    Long_t Nwgt, Nwgt_Negative;
    Double_t Swgt, Swgt_Negative;
    Double_t SQwgt, SQwgt_Negative;
    Long_t Nwgt_OverMax;
    Double_t WmaxTrue;
    Double_t Xsec, XsecErr;
    Double_t Xsec_Negative, Xsec_Negative_Err;
    Double_t Xsec_OverMax, Xsec_OverMax_Err;
   
    MCstat(std::ifstream &, bool debug=false);

    MCstat():
    Nevgen(0), Nwgt(0), Nwgt_Negative(0),
    Swgt(0), Swgt_Negative(0), 
    SQwgt(0), SQwgt_Negative(0),
    Nwgt_OverMax(0),
    WmaxTrue(0),
    Xsec(0), XsecErr(0),
    Xsec_Negative(0), Xsec_Negative_Err(0),
    Xsec_OverMax(0), Xsec_OverMax_Err(0)
    {};

    void SetEndofRun_mesmer();
    
    void Print() const;

    virtual ~MCstat(){};
    
    ClassDef(MCstat,3)
  };


  class Setup {

  public:
    MCpara MCpargen; 
    MCstat MCsums;

    Setup(const MCpara & p, const MCstat & s) : MCpargen(p), MCsums(s) {};
    Setup() : MCpargen(), MCsums() {};

    inline MCpara GetMCpara() const {return MCpargen;}
    virtual ~Setup(){};

    ClassDef(Setup,6)
  };

 class Particle {
  public:
    Short_t pdgId;
    Double_t px;
    Double_t py;
    Double_t pz;

    Double_t E() const;
    Double_t M() const;
    void Print() const;

    static void set_mass_e(const Double_t & me) {mass_e = me;}
    static void set_mass_mu(const Double_t & mm) {mass_mu = mm;}

    Particle(std::istringstream &);
    Particle():
    pdgId(0),px(0),py(0),pz(0)
    {};

    virtual ~Particle(){};
    
  private:
    static Double_t mass_e;
    static Double_t mass_mu;

    ClassDef(Particle,2)
  };

  class Event {

  public:
    Int_t RunNr;
    Long_t EventNr;
    Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO, wgt_NLO;
    static const uint NC = 11;
    std::array<Double_t, NC> coef;
    Particle muin;
    std::vector<Particle> fspart;

    Event():
    RunNr(0),EventNr(0),wgt_full(0),wgt_norun(0),wgt_lep(0),wgt_LO(0),wgt_NLO(0)
    {};

    bool Read(std::ifstream &, bool read_coef=true, bool debug=false);
    void Print() const;

    int GenerateEvent_mesmer(double* pmu);
    
    virtual ~Event(){};
    
    ClassDef(Event,7)
  };

}

#endif
