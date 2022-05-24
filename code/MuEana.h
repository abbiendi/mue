#ifndef MuEana_H
#define MuEana_H

///////////////////////////////////////////////
// Classes defining MuE analysis variables
//
// G.Abbiendi  4/Dec/2018 
///////////////////////////////////////////////

#include "Math/Vector4D.h"

namespace MuE {

  typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;

  // Analysis variables in the output tree
  
  class MuonOut {
  public:
    Double_t Emu; // muon energy
    Double_t thmu; // muon theta (in mrad)
    Double_t phmu; // muon phi (from -pi to +pi) 
    Double_t t13; // Mandelstam t (muon leg)
    Double_t x13; // Feynman x (muon leg)
    MuonOut():
    Emu(0),thmu(0),phmu(0),t13(0),x13(0)
      {};
    virtual ~MuonOut(){};
    ClassDef(MuonOut,1)
  };

  class Electron {
  public:
    Double_t Ee; // electron energy
    Double_t the; // electron theta (in mrad)
    Double_t phe; // electron phi (from -pi to +pi)
    Double_t t24; // Mandelstam t (electron leg)
    Double_t x24; // Feynman x (electron leg)
    Double_t tt_e; // t computed from electron angle with LO formulas
    Double_t xt_e; // x computed from electron angle with LO formulas
    Electron():
    Ee(0),the(0),phe(0),t24(0),x24(0),tt_e(0),xt_e(0)
      {};
    virtual ~Electron(){};
    ClassDef(Electron,1)
  };

  class MuEpair {
  public:
    Double_t deltaPhi; // acoplanarity (deltaPhi)
    Double_t openingAngle; // opening angle mu-e out in the Lab
    Double_t tripleProduct; // triple product btw normalized vectors i . mu x e
    MuEpair():
    deltaPhi(0),openingAngle(0),tripleProduct(0)
      {};
    virtual ~MuEpair(){};
    ClassDef(MuEpair,1)
  };

  class KineVars {
  public:
    Double_t t13; // Mandelstam t (muon leg)
    Double_t t24; // Mandelstam t (electron leg)
    Double_t x13; // Feynman x (muon leg)
    Double_t x24; // Feynman x (electron leg)
    Double_t tt_e; // t computed from electron angle with LO formulas
    Double_t xt_e; // x computed from electron angle with LO formulas

    Double_t Ee; // electron energy
    Double_t Emu; // muon energy
    Double_t the; // electron theta (in mrad)
    Double_t thmu; // muon theta (in mrad)
    Double_t phe; // electron phi (from -pi to +pi)
    Double_t phmu; // muon phi (from -pi to +pi) 

    Double_t deltaPhi; // acoplanarity (deltaPhi)
    Double_t openingAngle; // opening angle mu-e out in the Lab
    Double_t tripleProduct; // triple product btw normalized vectors i . mu x e

    void SetMuonOut(const MuonOut & m) {
      Emu = m.Emu;
      thmu = m.thmu;
      phmu = m.phmu;
      t13 = m.t13;
      x13 = m.x13;
    };
    void SetElectron(const Electron & e) {
      Ee = e.Ee;
      the = e.the;
      phe = e.phe;
      t24 = e.t24;
      x24 = e.x24;
      tt_e = e.tt_e;
      xt_e = e.xt_e;
    };
    void SetMuEpair(const MuEpair & d) {
      deltaPhi = d.deltaPhi;
      openingAngle = d.openingAngle;
      tripleProduct = d.tripleProduct;
    };

    KineVars():
      t13(0),t24(0),x13(0),x24(0),tt_e(0),xt_e(0),Ee(0),Emu(0),the(0),thmu(0),phe(0),phmu(0),
      deltaPhi(0),openingAngle(0),tripleProduct(0)
      {};
    
    KineVars(const MuonOut & m, const Electron & e, const MuEpair & d):
      t13(m.t13),t24(e.t24),x13(m.x13),x24(e.x24),tt_e(e.tt_e),xt_e(e.xt_e),
      Ee(e.Ee),Emu(m.Emu),the(e.the),thmu(m.thmu),phe(e.phe),phmu(m.phmu),
      deltaPhi(d.deltaPhi),openingAngle(d.openingAngle),tripleProduct(d.tripleProduct)
      {};

    virtual ~KineVars(){};
    ClassDef(KineVars,2)
  };

  class Photon {
  public:
    Double_t energy;    // photon energy in the Lab frame
    Double_t theta;     //   "    theta in the Lab frame (in mrad)
    Double_t phi;       //   "    phi in the Lab frame (in rad)
    Double_t energyCoM; // photon energy in the Centre-of-Mass frame
    
    Photon():
    energy(-1),theta(-1),phi(0),energyCoM(-1)
      {};
    
    virtual ~Photon(){};
    ClassDef(Photon,1)
  };

  class MuEana {

   public:
    Int_t RunNr;
    Long_t EventNr;
    Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO, wgt_NLO;   // event weights 
    PxPyPzEVector p_mu_in;
    PxPyPzEVector p_e_in;
    PxPyPzEVector p_mu_out;
    PxPyPzEVector p_e_out;
    std::vector<PxPyPzEVector> V_photons;
    KineVars genKin;   // kinematic variables at Generator-level for e and mu tracks
    KineVars detKin;   // kinematic variables at Detector-level for e and mu tracks
    Photon photon;     // photon kinematic variables at Gen-level
    
    MuEana():
    RunNr(0),EventNr(0),wgt_full(0),wgt_norun(0),wgt_lep(0),wgt_LO(0),wgt_NLO(0)
       {};
    
    virtual ~MuEana(){};
    ClassDef(MuEana,3)
  };
}

#endif
