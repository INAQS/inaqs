#include "properties.hpp"
#include "constants.hpp"
#include "bomd.hpp"
#include "qm_qchem.hpp"

BOMD::BOMD(size_t nqm, const int *qmid){
  qm = new QM_QChem(std::vector<int>(qmid, qmid + nqm), 0, 1);
  qm_grd.resize(nqm * 3);
  energy.resize(1);
};


template<typename T>
T BOMD::get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift)
{
  std::vector<double> qm_crd_v(qm_crd, qm_crd + 3*qm->nqm());
  std::vector<double> mm_crd_v(mm_crd, mm_crd + 3*nmm);
  std::vector<double> mm_chg_v(mm_chg, mm_chg + nmm);

  for(auto &x: qm_crd_v){x *= 10.0;} // nm->\AA
  for(auto &x: mm_crd_v){x *= 10.0;} // nm->\AA
  
  qm->update(qm_crd_v, mm_crd_v, mm_chg_v);

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

template double BOMD::get_gradient(const double* qm_crd, size_t nmm, const double* mm_crd, const double* mm_chg, double* f, double* fshift);
template float BOMD::get_gradient(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift);

template<typename T>
void BOMD::rescale_velocities(T* total_gradient, T* masses, T* velocities) {
  (void) total_gradient;
  (void) masses;
  (void) velocities;
};
template void BOMD::rescale_velocities(double* total_gradient, double* masses, double* velocities);
template void BOMD::rescale_velocities(float* total_gradient, float* masses, float* velocities);
