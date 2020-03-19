#ifndef GIFS_SH_BOMD_CORE_H
#define GIFS_SH_BOMD_CORE_H


#include <armadillo>
#include "qm_interface.hpp"


class BOMD
{
public:
  explicit BOMD(arma::uvec& atomicnumbers,
		arma::mat& qm_crd,
		arma::mat& mm_crd,
		arma::vec& mm_chg,
		arma::mat& qm_grd,
		arma::mat& mm_grd);
  virtual ~BOMD() {delete qm;}
  
  virtual double update_gradient(void);
  virtual bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift);

protected:
  // // fixed size
  // arma::cube qm_grd;
  // // flexible size
  // arma::cube mm_grd;
  // arma::vec energy;
    QMInterface* qm{nullptr};
    // fixed size
    arma::mat& qm_grd;
    // flexible size
    arma::mat& mm_grd;
    //
    arma::vec energy{};
};

class PrintBomd:
    public BOMD 
{
public:
  explicit PrintBomd(arma::uvec& atomicnumbers,
		arma::mat& qm_crd,
		arma::mat& mm_crd,
		arma::vec& mm_chg,
		arma::mat& qm_grd,
		arma::mat& mm_grd) : BOMD(atomicnumbers, qm_crd, mm_crd, mm_chg, qm_grd, mm_grd) {}
  
  double update_gradient(void) {
      qm->crd_qm.print("qm_crd: ");
      qm->crd_mm.print("mm_crd: ");
      qm_grd.fill(0.0);
      mm_grd.fill(0.0);
      return 0.0;
  }
  bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift) { 
      velocities.print("Velocities:");
      masses.print("Masses:");
      total_gradient.print("Total Gradient");
      return true;
  };
    

};


#endif
