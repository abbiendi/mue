#ifndef FastSim_H
#define FastSim_H

///////////////////////////////////////////////
// Fast Simulation of MuE scattering
//
// G.Abbiendi  4/Sep/2018 
///////////////////////////////////////////////
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/Vector2D.h"
#include "Math/Point3D.h"
#include "Math/Point2D.h"
#include "MuEtree.h"
#include "MuEana.h"
#include "Inputs.h"
#include <tuple>
#include "TMath.h"
#include <iostream>

namespace MuE {

typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;
typedef ROOT::Math::XYZPoint XYZPoint;
typedef ROOT::Math::XYZVector XYZVector;
typedef ROOT::Math::XYVector XYVector;

  class FastSim {
    
  public:
    FastSim(const MCpara & pargen, const FS_Input & parsim, bool _debug_=false);
    virtual ~FastSim(){};

    void Process(const Event & event);

    const PxPyPzEVector & Get_p_mu_in() const {return p_mu_in;}
    const PxPyPzEVector & Get_p_e_in() const {return p_e_in;}
    const PxPyPzEVector & Get_p_mu_out() const {return p_mu_out;}
    const PxPyPzEVector & Get_p_e_out() const {return p_e_out;}
    const std::vector<PxPyPzEVector> & Get_photons() const {return photons;}

    const KineVars & GetGenKin() const {return genKin;}
    const KineVars & GetDetKin() const {return detKin;}
    const Photon & GetPhoton() const {return photon;}

    void RandomNrSync();

  private:

    Double_t P_2bodies_CoM(Double_t Mass, Double_t mass_mu, Double_t mass_e) const;

    PxPyPzEVector Lorentz_ToCoM(const PxPyPzEVector & plab) const;
    PxPyPzEVector Lorentz_ToLab(const PxPyPzEVector & pcm) const;

    Double_t ThetaRMS(const PxPyPzEVector & p) const; 
    PxPyPzEVector Smear(const PxPyPzEVector & p, int istep=1) const; 
    PxPyPzEVector SmearX(const PxPyPzEVector & p) const; 
    PxPyPzEVector SmearPolar(const PxPyPzEVector & p) const; 
    PxPyPzEVector z_StationLength(const PxPyPzEVector & p) const; 
    PxPyPzEVector DetEffect(const PxPyPzEVector & p) const;
    KineVars LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
			  const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out);
    Photon LoadPhoton(const std::vector<PxPyPzEVector> & photons);
    MuonOut DetMuonOut(const PxPyPzEVector & p_mu_in, const PxPyPzEVector & p_mu_out) const;
    Electron DetElectron(const PxPyPzEVector & p_e_in, const PxPyPzEVector & p_e_out) const;
    MuEpair DetMuEpair(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
		       const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out) const;
        
    static const Double_t mm_PDG; // PDG muon mass 
    static const Double_t me_PDG; // PDG electron mass

    static const Double_t station_Length; // station length (1 m)
    static const Double_t detector_Size;  // transverse size of Si sensors (10 cm)
    
    const Double_t & mm; // muon mass
    const Double_t & me; // electron mass
    const Double_t & Qbeam; // beam muon charge
    const Double_t & Ebeam; // average beam energy
    const Double_t & EbeamRMS; // beam energy spread

    const Int_t & model; // model for detector smearing (gaussian resolution)
    const Int_t & MSopt; // options for multiple scattering (default/only Xplane/only polar)
    bool twosteps; // if true in model 0 apply only multiple scattering
    const Double_t & thickness; // material thickness (in X0) for model_=0
    const Double_t & intrinsic_resolution; // intrinsic resolution (in mrad) for model_=0
    const Double_t & zBias; // signed bias of station length (um)
    const Double_t & zSigma; // uncertainty of station length (um)
    bool zSigma_switch; // true if zSigma>0 
    bool debug;

    PxPyPzEVector p_system; // mu-e centre-of-mass system fourmomentum
    Double_t Minv; // event invariant mass = sqrt(s)

    PxPyPzEVector p_mu_in;
    PxPyPzEVector p_e_in;
    PxPyPzEVector p_mu_out;
    PxPyPzEVector p_e_out;
    std::vector<PxPyPzEVector> photons;

    Double_t thetaMax_ideal; // ideal max theta (geometric acceptance)
    Double_t thetaMax_cor;   // corrected max theta (geometric acceptance) for z bias

    KineVars genKin; // kinematic variables at Gen-level for e and mu track
    KineVars detKin; // kinematic variables at Detector-level for e and mu track
    Photon photon; // photon variables at Gen-level

  };
}

#endif
