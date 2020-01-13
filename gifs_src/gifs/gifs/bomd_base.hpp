#ifndef GIFS_SH_BOMD_BASE_H
#define GIFS_SH_BOMD_BASE_H

#include "properties.hpp"
#include "constants.hpp"

BOMD::BOMD(size_t nqm, std::vector<int>& qmid){
    qm = new QMInterface(nqm, qmid);
    qm_grd.resize(nqm * 3);
    energy.resize(1);
};


template<typename T>
T BOMD::get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift)
{
    qm->update(qm_crd, nmm, mm_crd, mm_chg);
    
    mm_grd.resize(nmm*3);

    PropMap props{};
    props.emplace(QMProperty::qmgradient, &qm_grd);
    props.emplace(QMProperty::mmgradient, &mm_grd);
    props.emplace(QMProperty::energies, &energy);

    qm->get_properties(props);

    // Unit conversion back qm->mm
    {
      int i = 0, j = 0;

      for (auto& fqm: qm_grd) {
        j = i++;
        f[j] = HARTREE_BOHR2MD*fqm;
        fshift[j] = f[j];
      }

      for (auto& fmm: mm_grd) {
        j = i++;
        f[j] = HARTREE_BOHR2MD*fmm;
        fshift[j] = f[j];
      }
    }
    // Return in "MD" units of KJ/mole 
    return HARTREE2KJ * AVOGADRO * energy[0];
};


template<typename T>
void BOMD::rescale_velocities(T* total_gradient, T* masses, T* velocities) {
  
};

#endif
