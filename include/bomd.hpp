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
    //
    virtual double update_gradient();
    //
    virtual ~BOMD() {delete qm;}
    //
    template<typename T> // templates cannot be virtual
    void rescale_velocities(T* total_gradient, T* masses, T* velocities);

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
