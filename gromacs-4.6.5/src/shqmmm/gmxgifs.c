#include <stdio.h>
#include "gmxgifs.h"
#include "gifs.hpp"

void gifs_scale_velocities(float *v, float *f, float* invmass) {
    printf("WE DO RESCALING!\n");
    gifs_rescale_velocities(f, invmass, v);
};

float gifs_print_coords(int nqm, const int* qm_atomids, const float* qm_crd) {
   int i;
  /*
  create_qm_interface(nqm, qm_atomids);
  return gifs_get_forces(qm_crd, nmm, mm_crd,  mm_chg, f, fshift);
  */
  printf("QM Coords: nqm = %d\n", nqm);
  for(i=0; i<nqm; ++i) {
    printf("%d %12.8f %12.8f %12.8f\n", qm_atomids[i], qm_crd[i*3]*10.0, qm_crd[i*3+1]*10.0, qm_crd[i*3+2]*10.0);
  }
  printf("DONE\n");
  return 0.0;
};

void gifs_update_global_idx(int* indexQM, int* indexMM) {
  gifs_update_global_index(indexQM, indexMM);
};

float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
			size_t nmm, const float* mm_crd, const float* mm_chg,
			float* f, float* fshift){
   int i;
   create_qm_interface("./gifs_config.ini", nqm, qm_atomids);
  return gifs_get_forces(qm_crd, NULL, nmm, mm_crd,  mm_chg, f, f + nqm*3);
};
