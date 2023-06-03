#ifndef GMX_GIFS_INTERFACE_H
#define GMX_GIFS_INTERFACE_H
#include <stddef.h>
#include <stdbool.h>
#include "typedefs.h"

void gifs_update_global_idx(int* indexQM, int* indexMM);
void gifs_scale_velocities(real energy, rvec *v, rvec *f, real* invmass);
real gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const rvec* qm_crd,
			size_t nmm, const rvec* mm_crd, const real* mm_chg,
                       rvec* f, rvec* fshift);

bool inaqs_init(const char * configfile, int mdStepStart, real timeStep, size_t nqm, const int * qm_atomids);

#endif
