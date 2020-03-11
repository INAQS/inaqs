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
  virtual void rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::vec &total_gradient, double e_drift);

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


#endif
