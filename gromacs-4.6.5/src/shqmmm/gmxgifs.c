#include <stdio.h>
#include "gmxgifs.h"
#include "gifs.hpp"

void gifs_scale_velocities(float* v, float* f, float* invmass) {}
float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
			size_t nmm, const float* mm_crd, const float* mm_chg,
			float* f, float* fshift){
  create_qm_interface(nqm, qm_atomids);
  return gifs_get_forces(qm_crd, nmm, mm_crd,  mm_chg, f, fshift);
}
