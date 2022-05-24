#include <cmath>
#include <iostream>
#include "TRandom.h"
#include "ElasticState.h"
#include "ResolutionModels.h"
#include "FastSim.h"

using namespace MuE;
using namespace std;
using namespace ROOT::Math;

// from PDG Book 2018
const Double_t FastSim::mm_PDG = 105.6583745 *0.001;
const Double_t FastSim::me_PDG = 0.5109989461 *0.001;

// detector geometry
const Double_t FastSim::station_Length = 1.; // 1m
const Double_t FastSim::detector_Size = 0.1; // 10cm

FastSim::FastSim(const MuE::MCpara & pargen, const MuE::FS_Input & parsim, bool _debug_):
  mm(pargen.mass_mu), me(pargen.mass_e), Qbeam(pargen.charge_mu), Ebeam(pargen.Ebeam), EbeamRMS(pargen.EbeamRMS),
  model(parsim.model), MSopt(parsim.MSopt), 
  twosteps(parsim.twosteps), thickness(parsim.thickness), intrinsic_resolution(parsim.resolution), 
  zBias(parsim.zBias), zSigma(parsim.zSigma), debug(_debug_), Minv(0)
  
{
  if (std::abs(mm - mm_PDG)/mm_PDG > 1e-7) {
    cout<<"\n"<< "***WARNING: muon mass = "<<mm<<" is different from the PDG mass: "<<mm_PDG<<endl;
  }
  if (std::abs(me - me_PDG)/me_PDG > 1e-7) {
    cout<<"\n"<< "***WARNING: electron mass = "<<me<<" is different from the PDG mass: "<<me_PDG<<endl;
  }

  if (model == 0) {
    cout<<"\n Simple detector model, DetModel = "<<model<<", option = "<<MSopt<<", twosteps = "<<twosteps <<endl;
    cout<<" material thickness = "<<thickness<<" X0, intrinsic angular resolution = "<<intrinsic_resolution<<" mrad "<<endl;
  }
  else if (model == 1) {
    cout<<"\n Antonios detector model, 3 parameters"<< endl;
    twosteps = false;
  }
  else {
    cout<<"\n"<<"*** ERROR : FastSim, undefined detector model with DetModel = "<<model<<endl; 
    std::exit(999);
  }
  
  thetaMax_ideal = detector_Size/2 /station_Length;
  thetaMax_cor = thetaMax_ideal*(1.-zBias*1e-6/station_Length);

  if (zSigma > 1e-7) {
    zSigma_switch = true;   
    cout <<" ===> setting the zSigma_switch to true"<<endl;
  }
  else {
    zSigma_switch = false;
    cout <<" ===> setting the zSigma_switch to false"<<endl;
  }

  cout<<"\n Station Z length : bias = "<<zBias<<" (um), sigma = "<<zSigma<<" (um) "<<endl;
  cout <<" thetaMax_ideal = "<<thetaMax_ideal*1e3 <<" mrad, thetaMax_cor = "<<thetaMax_cor*1e3 <<" mrad"<<endl;
}


// process an event
//
void FastSim::Process(const MuE::Event & event) {
  
  if (debug) cout<<"\n Process:  Run = "<<event.RunNr << " , Event = "<< event.EventNr << endl;
 
  p_mu_in.SetPxPyPzE(event.muin.px, event.muin.py, event.muin.pz, event.muin.E());
  p_e_in.SetPxPyPzE(0, 0, 0, me);

  p_system = p_mu_in + p_e_in;
  Double_t s = mm*mm + me*me + 2*me*p_mu_in.E();
  Minv = sqrt(s);

  //  Double_t pcm = P_2bodies_CoM(Minv, mm, me);
  if (debug) cout<<"\n"<<"Incoming muon energy = "<<std::fixed <<setprecision(3)<< p_mu_in.E() <<" GeV"<<endl;
  //      << " GeV, s = "<<s<<" GeV^2, sqrt(s) = "<<Minv<<" GeV, pcm = "<<pcm<<" GeV"<<endl;
  
  bool found_mu = false;
  bool found_e = false;
  UInt_t nphot = 0;
  photons.clear();

  for (const auto & p : event.fspart) {
    if (!found_mu && ((p.pdgId) == -Qbeam*13)) {
      found_mu = true;
      p_mu_out.SetPxPyPzE(p.px, p.py, p.pz, p.E());
    }
    else if (!found_e && (p.pdgId == 11)) {
      found_e = true;
      p_e_out.SetPxPyPzE(p.px, p.py, p.pz, p.E());
    }
    else if (nphot<2 && p.pdgId == 22) {
      nphot++;
      photons.push_back(PxPyPzEVector(p.px, p.py, p.pz, p.E()));
    }
    else {
      cout<<"\n*** ERROR: unexpected particle in Event "<<event.EventNr <<endl;
      cout<<"\t incoming muon: [id:"<<event.muin.pdgId<<"] ("<<event.muin.px<<", "<<event.muin.py<<", "<<event.muin.pz<<")"<<endl;
      cout<<"\t final-state particles: "<<endl;
      for (const auto & fsp : event.fspart)
	cout<<"\t [id:"<<fsp.pdgId<<"] ("<<fsp.px<<", "<<fsp.py<<", "<<fsp.pz<<")"<<endl;
      std::exit(222);
    }
  }
  
  // Generator-level kinematics
  genKin = LoadKineVars(p_mu_in, p_e_in, p_mu_out, p_e_out);

  // smearing (multiple scattering on target and optionally intrinsic resolution on the detector)
  PxPyPzEVector p_mu_out_sim, p_e_out_sim; 
  if (MSopt ==0) {
    p_mu_out_sim = Smear(p_mu_out);
    p_e_out_sim = Smear(p_e_out);
  }
  else if (MSopt ==1) {
    p_mu_out_sim = SmearX(p_mu_out);
    p_e_out_sim = SmearX(p_e_out);
  }
  else if (MSopt ==2) {
    p_mu_out_sim = SmearPolar(p_mu_out);
    p_e_out_sim = SmearPolar(p_e_out);
  }
  else cout<<"\n"<<"*** ERROR : FastSim, undefined detector MS option = "<<MSopt<<endl;
  
  // no smearing yet for incoming muon
  PxPyPzEVector p_mu_in_det = p_mu_in;
  // Detector effects: geometric acceptance and intrinsic resolution
  PxPyPzEVector p_mu_out_det = DetEffect(p_mu_out_sim);
  PxPyPzEVector p_e_out_det = DetEffect(p_e_out_sim); 
  
  // final printouts
  if (debug) cout<<"\n=== GEN Level === "<<endl;
  if (debug) cout<<"p_mu_in = "<<p_mu_in.Px()<<", "<<p_mu_in.Py()<<", "<<p_mu_in.Pz()<< endl;
  if (debug) cout<<"p_mu_out = "<<p_mu_out.Px()<<", "<<p_mu_out.Py()<<", "<<p_mu_out.Pz()<< endl;
  if (debug) cout<<"p_e_out = "<<p_e_out.Px()<<", "<<p_e_out.Py()<<", "<<p_e_out.Pz()<< endl;
  if (debug) cout<<"theta[mrad] (in)= "<<1e3*p_mu_in.Theta()
      <<", (mu)= "<<1e3*p_mu_out.Theta()<<", (el)= "<<1e3*p_e_out.Theta()<<endl;
  // load final detector-level quantities
  if (debug) cout<<"\n=== SIM-RECO === "<<endl;
  if (debug) cout<<"p_mu_in_det = "<<p_mu_in_det.Px()<<", "<<p_mu_in_det.Py()<<", "<<p_mu_in_det.Pz()<< endl;
  if (debug) cout<<"p_mu_out_det = "<<p_mu_out_det.Px()<<", "<<p_mu_out_det.Py()<<", "<<p_mu_out_det.Pz()<< endl;
  if (debug) cout<<"p_e_out_det = "<<p_e_out_det.Px()<<", "<<p_e_out_det.Py()<<", "<<p_e_out_det.Pz()<< endl;
  if (debug) cout<<"theta[mrad] (in)= "<<1e3*p_mu_in_det.Theta()
      <<", (mu)= "<<1e3*p_mu_out_det.Theta()<<", (el)= "<<1e3*p_e_out_det.Theta()<<endl;

  detKin = LoadKineVars(p_mu_in_det, p_e_in, p_mu_out_det, p_e_out_det);
  photon = LoadPhoton(photons); // well, this is still gen-level with no acceptance cut!

}


KineVars FastSim::LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
			       const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out) 
{
  MuonOut m = p_mu_out.E() > 0 ? DetMuonOut(p_mu_in, p_mu_out) : MuonOut();
  Electron e = p_e_out.E() > 0 ? DetElectron(p_e_in, p_e_out) : Electron();
  MuEpair d = (p_mu_out.E() > 0 && p_e_out.E() > 0) ? DetMuEpair(p_mu_in, p_e_in, p_mu_out, p_e_out) : MuEpair();

  return KineVars(m,e,d);
}


MuonOut FastSim::DetMuonOut(const PxPyPzEVector & p_mu_in, const PxPyPzEVector & p_mu_out) const
{
  MuonOut m;
  m.Emu = p_mu_out.E();
  m.thmu = 1e3* p_mu_out.Theta();
  m.phmu = p_mu_out.Phi();
  
  PxPyPzEVector q13 = p_mu_in - p_mu_out;
  m.t13 = q13.M2();
  ElasticState emu(Ebeam,mm,me);
  m.x13 = emu.X(m.t13);

  return m;
}


Electron FastSim::DetElectron(const PxPyPzEVector & p_e_in, const PxPyPzEVector & p_e_out) const
{
  Electron e;
  e.Ee = p_e_out.E();
  e.the = 1e3* p_e_out.Theta();
  e.phe = p_e_out.Phi();
  
  // Note: here Ebeam is the average beam energy, so tt_e and xt_e are defined under this assumption
  MuE::ElasticState emu(Ebeam,mm,me, e.the);
  e.tt_e = emu.GetT();
  e.xt_e = emu.GetX();
  
  PxPyPzEVector q24 = p_e_in - p_e_out;
  // instead these are exact
  e.t24 = q24.M2();
  e.x24 = emu.X(e.t24);

  return e;
}


MuEpair FastSim::DetMuEpair(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
			    const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out) const
{  
  MuEpair d;

  // acoplanarity = deltaPhi as long as the incoming muon is collinear with z-axis
  Double_t deltaPhi = p_e_out.Phi() - p_mu_out.Phi();
  if (deltaPhi<0.) deltaPhi = deltaPhi + 2.*TMath::Pi();
  deltaPhi = deltaPhi - TMath::Pi();
  d.deltaPhi = deltaPhi;
  
  XYZVector nvec_mu_in = p_mu_in.Vect().Unit(); 
  XYZVector nvec_mu_out = p_mu_out.Vect().Unit();
  XYZVector nvec_e_out = p_e_out.Vect().Unit();

  Double_t dotProduct = nvec_mu_out.Dot(nvec_e_out);
  d.openingAngle = std::abs(dotProduct)<1. ? 1000.*acos(dotProduct) : 0.;

  XYZVector crossProduct = nvec_mu_out.Cross(nvec_e_out);
  d.tripleProduct = nvec_mu_in.Dot(crossProduct);

  return d;
}


// RMS Theta smearing due to Multiple scattering distribution and intrinsic resolution
//
Double_t FastSim::ThetaRMS(const PxPyPzEVector & k) const
{
  Double_t pmom = k.P();
  
  // resolution: gaussian sigma
  Double_t thrms(0); 
  
  if (model == 0) {
    Double_t msc = thickness > 0 ? MuE::Target_thrms(pmom, thickness) : 0; // in rad
    // N.B. intrinsic_resolution is in mrad
    thrms = twosteps ? msc : sqrt(msc*msc + 1e-6*intrinsic_resolution*intrinsic_resolution);
  }
  else if (model == 1) {
    thrms = MuE::Antonio_thrms(pmom); // in mrad
  }
  else {
    cout << "\n" << "***ERROR: Undefined smearing model = "<<model << endl;
    exit(999);
  }

  return thrms;
}

// apply resolution smearing to particle momentum
//
PxPyPzEVector FastSim::Smear(const PxPyPzEVector & k, int istep) const
{
  
  Double_t thrms = istep==1 ? ThetaRMS(k) : intrinsic_resolution;
  
  Double_t smearx = gRandom->Gaus(0., thrms);
  Double_t smeary = gRandom->Gaus(0., thrms);
  
  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??
  
  // apply smearing
  anglex += smearx;
  angley += smeary;
  
  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);
  
  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz*dxdz + dydz*dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;
  
  PxPyPzEVector psmeared(skx, sky, skz, k.E());
  
  return psmeared;
}


PxPyPzEVector FastSim::DetEffect(const PxPyPzEVector & p_sim) const
{  
  PxPyPzEVector p = p_sim;
  
  // apply z-scale and uncertainty of station length
  if (zBias != 0 || zSigma != 0) p = z_StationLength(p_sim); 
  
  // check geometric acceptance
  if (p.Theta() > thetaMax_cor) {
    //cout<<"theta = "<<p.Theta() <<" > "<<thetaMax_cor<<":  twosteps = "<< twosteps <<endl;
    if (twosteps) {
      //      cout<<"\t calling two randoms "<<endl;
      gRandom->Gaus(0., 1.); // dummy calls to preserve the random chain synchronization
      gRandom->Gaus(0., 1.); 
    }
    return PxPyPzEVector(0,0,0,0);
  }
  
  // apply intrinsic detector resolution for two-steps model 0
  if (twosteps) {
    p = Smear(p, 2);
    //   cout <<"intrinsic smearing " << endl;
  }

  return p;
}


PxPyPzEVector FastSim::z_StationLength(const PxPyPzEVector & k) const
{
  Polar3DVector p(k.P(), k.Theta(), k.Phi()); 
  double the = p.Theta();

  double relbias = (zBias*1e-6) /station_Length; // zBias is in um 
  double thebias = -the*relbias;
  double thesm = the + thebias;

  double thesigma = the * (zSigma*1e-6) /station_Length;
  
  if (zSigma_switch) {
    thesm = thesm +  gRandom->Gaus(0., thesigma);
  }
  
  double phism = p.Phi();

  // note that theta is defined positive (0-pi)
  // going negative means changing the azimuth too and phi is defined in (-pi,+pi)
  if (thesm < 0) {
    thesm = -thesm;
    phism = phism > 0 ? phism - TMath::Pi() : phism + TMath::Pi();
    cout<<"*** z_StationLength: thesm < 0  = "<< thesm*1e3 <<" mrad" <<endl;
  }

  p.SetTheta(thesm);
  p.SetPhi(phism);

  return PxPyPzEVector(p.X(), p.Y(), p.Z(), k.E());
}

/*
void FastSim::LoadPhoton(const MuE::Event & event, MuE::Photon & photon) {
  // by now at most one photon
  auto n_photons = event.photons.size();
  
  if (n_photons >0) { 
    PxPyPzEVector p_gamma_Lab = {event.photons[0].px, 
				 event.photons[0].py,
				 event.photons[0].pz,
			         event.photons[0].E};
    PxPyPzEVector p_gamma_CoM = Lorentz_ToCoM(p_gamma_Lab);
    
    photon.energy    = p_gamma_Lab.E();
    photon.theta     = p_gamma_Lab.Theta() *1e3;
    photon.phi       = p_gamma_Lab.Phi();
    photon.energyCoM = p_gamma_CoM.E();  
  }

  else {
    photon.energyCoM = -1;
    photon.energy    = -1;
    photon.theta     = -1;
    photon.phi       =  0;
  }
}
*/
Photon FastSim::LoadPhoton(const std::vector<PxPyPzEVector> & photons) {
  // NLO: at most one photon
  auto n_photons = photons.size();
  
  MuE::Photon photon;
  if (n_photons >0) { 
    PxPyPzEVector p_gamma_Lab = photons[0];
    PxPyPzEVector p_gamma_CoM = Lorentz_ToCoM(p_gamma_Lab);
    
    photon.energy    = p_gamma_Lab.E();
    photon.theta     = p_gamma_Lab.Theta() *1e3;
    photon.phi       = p_gamma_Lab.Phi();
    photon.energyCoM = p_gamma_CoM.E();  
  }
  return photon;
}

// momentum for 2-body kinematics in the centre-of-mass system 
//
Double_t FastSim::P_2bodies_CoM (Double_t M, Double_t mm, Double_t me) const
{
  Double_t msum = std::abs(mm+me);
  Double_t mdif = std::abs(mm-me); 
  return 0.5/M*sqrt((M+mdif)*(M-mdif)*(M+msum)*(M-msum));
}

// Lorentz transformation to the centre-of-mass system
//
PxPyPzEVector FastSim::Lorentz_ToCoM(const PxPyPzEVector & pLab) const
{ 
  if (p_system.E() == Minv) return pLab;
  
  else {
    Double_t ecm = (pLab.E()*p_system.E()
		    -pLab.Px()*p_system.Px()-pLab.Py()*p_system.Py()-pLab.Pz()*p_system.Pz()) /Minv;
    
    Double_t fn = (ecm+pLab.E()) / (p_system.E()+Minv);
    
    return PxPyPzEVector(pLab.Px() - fn * p_system.Px(),
			  pLab.Py() - fn * p_system.Py(),
			  pLab.Pz() - fn * p_system.Pz(),
			  ecm);
  }
}

// Lorentz transformation to the Laboratory system
//
PxPyPzEVector FastSim::Lorentz_ToLab(const PxPyPzEVector & pCoM) const
{
  if (p_system.E() == Minv) return pCoM;
  
  else {
    Double_t elab = (pCoM.Px()*p_system.Px()+pCoM.Py()*p_system.Py()
		     +pCoM.Pz()*p_system.Pz()+pCoM.E()*p_system.E()) /Minv;
    
    Double_t fn   = (elab+pCoM.E()) / (p_system.E()+Minv);
    
    return PxPyPzEVector(pCoM.Px() + fn * p_system.Px(),
			  pCoM.Py() + fn * p_system.Py(),
			  pCoM.Pz() + fn * p_system.Pz(),
			  elab);
  }
}


// apply resolution smearing to particle momentum (SMEAR only on the XZ plane)
//
PxPyPzEVector FastSim::SmearX(const PxPyPzEVector & k) const
{
  Double_t thrms = ThetaRMS(k);
  
  Double_t smearx = gRandom->Gaus(0., thrms);
  gRandom->Gaus(0., 1.); // dummy call to preserve the random chain synchronization
  Double_t smeary = 0;
  
  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??
  
  // apply smearing
  anglex += smearx;
  angley += smeary;
  
  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);
  
  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz*dxdz + dydz*dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;
  
  PxPyPzEVector psmeared(skx, sky, skz, k.E());
  
  return psmeared;
}

// apply resolution smearing to particle momentum
// (theta smearing of rms*sqrt(2) in the plane defined by the vector and the z-axis)
//
PxPyPzEVector FastSim::SmearPolar(const PxPyPzEVector & k) const
{
  Polar3DVector p(k.P(), k.Theta(), k.Phi()); 

  // assumed smearing is sqrt(2) * smearing on a plane
  Double_t thrms = ThetaRMS(k);  
  Double_t smearth = sqrt(2) * gRandom->Gaus(0., thrms);
  gRandom->Gaus(0., 1.); // dummy call to preserve the random chain synchronization
  Double_t thetasm = p.Theta() + smearth;
  Double_t phism = p.Phi();
  // note that theta is defined positive (0-pi)
  // going negative means changing the azimuth too and phi is defined in (-pi,+pi)
  if (thetasm < 0) {
    thetasm = -thetasm;
    phism = phism > 0 ? phism - TMath::Pi() : phism + TMath::Pi();
  }

  p.SetTheta(thetasm);
  p.SetPhi(phism);

  return PxPyPzEVector(p.X(), p.Y(), p.Z(), k.E());
}


// synchronize the random number chain to account for events with negligible weight
//  (skipped in the main event loop)
//
void FastSim::RandomNrSync() {
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing muon theta X
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing muon theta Y
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing electron theta X
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing electron theta Y
}
