#ifndef GMX_GIFS_INTERFACE_H
#define GMX_GIFS_INTERFACE_H

void gifs_scale_velocities(float* v, float* f, float* invmass);
float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
			size_t nmm, const float* mm_crd, const float* mm_chg,
                       float* f, float* fshift);
#endif
