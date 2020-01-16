#ifndef GIFS_MAIN_H
#define GIFS_MAIN_H

#ifdef __cplusplus
extern "C" {
#endif
 
  void create_qm_interface(size_t nqm, const int* qm_atomids);
  float gifs_get_forces(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift);

#ifdef __cplusplus 
}
#endif

#endif
