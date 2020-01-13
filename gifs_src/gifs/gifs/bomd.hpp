#ifndef GIFS_SH_BOMD_CORE_H
#define GIFS_SH_BOMD_CORE_H

#include <vector>

class BOMD
{
public:

    explicit BOMD(size_t nqm, std::vector<int>& qmid);

    ~BOMD() {delete qm;}

    template<typename T>
    T get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift);
    //
    template<typename T>
    inline T rescale_velocities(T* total_gradient, T* masses, T* velocities);

protected:
    QMInterface* qm{nullptr};
    // fixed size
    std::vector<double> qm_grd{};
    // flexible size
    std::vector<double> mm_grd{};
    std::vector<double> energy{};
};



#endif
