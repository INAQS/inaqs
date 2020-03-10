#include "properties.hpp"
#include "constants.hpp"
#include "bomd.hpp"
#include "qm_qchem.hpp"
#include <armadillo>


BOMD::BOMD(arma::uvec& atomicnumbers,
           arma::mat& qm_crd, 
           arma::mat& mm_crd, 
           arma::vec& mm_chg, 
           arma::mat& qm_grd,
           arma::mat& mm_grd) :
    qm_grd{qm_grd}, mm_grd{mm_grd}, energy(1)
{
  qm = new QM_QChem(atomicnumbers, qm_crd, mm_crd, mm_chg, 0, 1, 0);
};


double 
BOMD::update_gradient()
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


void BOMD::rescale_velocities(arma::vec &velocities, arma::vec &masses, arma::vec &total_gradient, double e_drift) {
  (void) total_gradient;
  (void) masses;
  (void) velocities;
  (void) e_drift;
};


template void BOMD::rescale_velocities(double* total_gradient, double* masses, double* velocities);
template void BOMD::rescale_velocities(float* total_gradient, float* masses, float* velocities);
