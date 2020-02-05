#include "properties.hpp"
#include "constants.hpp"
#include "bomd.hpp"
#include "qm_qchem.hpp"
#include <armadillo>


// BOMD::BOMD(size_t nqm, const int *qmid){
//   qm = new QM_QChem(std::vector<int>(qmid, qmid + nqm), 0, 1, 0);
//   qm_grd.resize(3, nqm, 1);
//   energy.resize(1);
BOMD::BOMD(std::vector<int.& qmids,
           std::vector<double>& qm_crd, 
           std::vector<double>& mm_crd, 
           std::vector<double>& mm_chg, 
           std::vector<double>& in_qm_grd,
           std::vector<double>& in_mm_grd) :
    qm_grd{in_qm_grd}, mm_grd{in_mm_grd}, energy{}
{
  qm = new QM_QChem(qmids, qm_crd, mm_crd, mm_chg, 0, 1);
>>>>>>> b87ce38... updated interface
};


double BOMD::get_gradient()
{
    qm->update();

    PropMap props{};
    props.emplace(QMProperty::qmgradient, &qm_grd);
    props.emplace(QMProperty::mmgradient, &mm_grd);
    props.emplace(QMProperty::energies, &energy);
    //
    qm->get_properties(props);
    //
    return energy[0];
};


template<typename T>
void BOMD::rescale_velocities(T* total_gradient, T* masses, T* velocities) {
  (void) total_gradient;
  (void) masses;
  (void) velocities;
};


template void BOMD::rescale_velocities(double* total_gradient, double* masses, double* velocities);
template void BOMD::rescale_velocities(float* total_gradient, float* masses, float* velocities);
