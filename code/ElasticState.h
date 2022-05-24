#ifndef ElasticState_H
#define ElasticState_H

///////////////////////////////////////////////
// Kinematics for MuE scattering
//
// G.Abbiendi  4/Sep/2018 
///////////////////////////////////////////////

namespace MuE {
  
  class ElasticState {
    
  public:
    ElasticState(double energy_mu_in, double mass_mu, double mass_e, double theta_e_out);
    ElasticState(double energy_mu_in, double mass_mu, double mass_e);
    ElasticState(){};    
    virtual ~ElasticState(){};
    
    double GetTheta_e() const {return the_;}
    double GetTheta_mu() const {return thmu_;}
    double GetEnergy_e() const {return Ee_;}
    double GetEnergy_mu() const {return Emu_;}
    double GetT() const {return t_;}
    double GetX() const {return x_;}

    //double X(const double & T) const;
    inline double X(const double & t) const {
      return (t + sqrt(t*t - 4*mm_*mm_*t))/(2*mm_*mm_);
    }

    //double tx(const double & X);
    // calculate t given x
    inline double tx(const double & x) const {
      return -mm_*x/(1.-x)*mm_*x;
    }
    
    //double Eet(const double & T) const;
    // calculate Ee (outgoing electron energy) given t
    inline double Eet(const double & t) const {
      return me_ - t/2./me_;
    }

    //double theEe(const double & energy_e_out) const;
    inline double theEe(const double & Ee) const {
      double beta = sqrt(1. - (mm_/E_)*(mm_/E_)) / (1. + me_/E_);
      return acos(sqrt((Ee-me_)/(Ee+me_)) /beta);
    }
    
  private:
    double Ee(const double & theta_e_out) const;
    double Emu(const double & energy_e_out) const;
    double Thmu(const double & energy_mu_out) const;
    double T(const double & energy_e_out) const;
    
    double E_, mm_, me_;
    double the_, thmu_, Ee_, Emu_, t_, x_;
  };
}

#endif
