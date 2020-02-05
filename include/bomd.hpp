#ifndef GIFS_SH_BOMD_CORE_H
#define GIFS_SH_BOMD_CORE_H


#include <vector>
#include <armadillo>
#include "qm_interface.hpp"


class BOMD
{
public:
    explicit BOMD(std::vector<int>& qmids,
                  std::vector<double>& qm_crd, 
                  std::vector<double>& mm_crd, 
                  std::vector<double>& mm_chg, 
                  std::vector<double>& qm_grd,
                  std::vector<double>& mm_grd);
    //
    double get_gradient();
    //
    virtual ~BOMD() {delete qm;}
    //
    template<typename T>
    void rescale_velocities(T* total_gradient, T* masses, T* velocities);

protected:
  // // fixed size
  // arma::cube qm_grd;
  // // flexible size
  // arma::cube mm_grd;
  // arma::vec energy;
    QMInterface* qm{nullptr};
    // fixed size
    std::vector<double>& qm_grd;
    // flexible size
    std::vector<double>& mm_grd;
    std::vector<double> energy{};
};


#endif
