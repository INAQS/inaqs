#ifndef GIFS_MAIN_H
#define GIFS_MAIN_H

#ifdef __cplusplus
#include <memory>
#include <armadillo>
#include "qm_interface.hpp"

extern "C" {
#endif
  void create_qm_interface(const char* file, int mdStepStart, double classicalTimeStep, int nqm, const int* qm_atomids);
  float gifs_get_forces_float(const float* qm_crd, const size_t* link_ids, size_t nmm, const float* mm_crd, const float* mm_chg, float* f_qm, float* f_mm);
  double gifs_get_forces_double(const double* qm_crd, const size_t* link_ids, size_t nmm, const double* mm_crd, const double* mm_chg, double* f_qm, double* f_mm);
  void gifs_update_global_index(int* indexQM, int* indexMM);
  void gifs_rescale_velocities_float(float energy, float* total_gradient, float* masses, float* velocities);
  void gifs_rescale_velocities_double(double energy, double* total_gradient, double* masses, double* velocities);
  bool gifs_interface_is_ready(void);
#ifdef __cplusplus 
}

std::shared_ptr<QMInterface> gifs_QMInterface(void);
void gifs_update_coords(const arma::mat & R);
#endif

#endif
