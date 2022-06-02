#include <stdio.h>
#include "gmxgifs.h"
#include "gifs.hpp"
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>

void gifs_scale_velocities(float energy, float *v, float *f, float* invmass) {
  gifs_rescale_velocities(energy, f, invmass, v);
};

void gifs_update_global_idx(int* indexQM, int* indexMM) {
  gifs_update_global_index(indexQM, indexMM);
};

float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
			size_t nmm, const float* mm_crd, const float* mm_chg,
			float* f, float* fshift){

  if (!gifs_interface_is_ready()){
    fprintf(stderr, "[INAQS] WARNING: attempt to call %s before init().\n", __func__);
  }

  return gifs_get_forces(qm_crd, NULL, nmm, mm_crd,  mm_chg, f, f + nqm*3);
};

bool inaqs_init(char * inaqsConfigFile, float timeStep, size_t nqm, const int * qm_atomids){
  (void) inaqsConfigFile;  // FIXME: need to actually use
  (void) timeStep;         // FIXME: need to actually use

  printf("NQM = %ld:", nqm);
  for (size_t i = 0; i < nqm; i++){
    printf(" %d", qm_atomids[i]);
  }
  printf("\n");
  
  char * config_names[] = {
    "inaqs_config.ini", // default name is first
    "gifs_config.ini",  // legacy inputs later
    NULL                // NULL-terminate
  };
  
  for (char ** fname = config_names; *fname; fname++){
    if (0 == access(*fname, F_OK)){
      if (fname != config_names){
        fprintf(stderr,
                "DEPRECATION WARNING: you are not using the standard config file"
                ", '%s'. Please update your inputs.\n", config_names[0]);
      }
      create_qm_interface(*fname, nqm, qm_atomids);
      return true;
    }
  }
    
  fprintf(stderr, "Unable to locate INAQS config file, '%s', terminating!\n",
          config_names[0]);
  exit(-1);
}
