#ifndef GIFS_MAIN_H
#define GIFS_MAIN_H

#ifdef __cplusplus
#include <memory>
#include "qm_interface.hpp"

extern "C" {
#endif
 
//void create_qm_interface(size_t nqm, const int* qm_atomids);
//float gifs_get_forces(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift);
//float gifs_get_forces(float* qm_crd, int nmm, float* mm_crd, float* mm_chg, float* f, float* fshift);
//void gifs_update_global_index(int* indexQM, int* indexMM);
//void gifs_rescale_velocities(float* total_gradient, float* masses, float* velocities);
void create_qm_interface(const char* file, int nqm, const int* qm_atomids);
float gifs_get_forces(const float* qm_crd, const size_t* link_ids, size_t nmm, const float* mm_crd, const float* mm_chg, float* f_qm, float* f_mm);
void gifs_update_global_index(int* indexQM, int* indexMM);
  void gifs_rescale_velocities(float energy, float* total_gradient, float* masses, float* velocities);

#ifdef __cplusplus 
}

std::shared_ptr<QMInterface> gifs_QMInterface(void);
#endif

#endif
