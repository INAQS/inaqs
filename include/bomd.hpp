#ifndef GIFS_SH_BOMD_CORE_H
#define GIFS_SH_BOMD_CORE_H

#include <vector>
#include <armadillo>
#include "qm_interface.hpp"

class BOMD
{
public:

    explicit BOMD(size_t nqm, const int *qmid);

    virtual ~BOMD() {delete qm;}

    template<typename T>
    T get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift);
    //
    template<typename T>
    void rescale_velocities(T* total_gradient, T* masses, T* velocities);

protected:
  QMInterface* qm{nullptr};
  // fixed size
  arma::cube qm_grd;
  // flexible size
  arma::cube mm_grd;
  arma::cube energy;
};



#endif
