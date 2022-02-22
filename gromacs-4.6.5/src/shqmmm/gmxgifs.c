#include <stdio.h>
#include "gmxgifs.h"
#include "gifs.hpp"
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>

void gifs_scale_velocities(float energy, float *v, float *f, float* invmass) {
  gifs_rescale_velocities(energy, f, invmass, v);
};

float gifs_print_coords(int nqm, const int* qm_atomids, const float* qm_crd) {
   int i;
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
  static bool needinit = true;

  if (needinit){
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
        needinit = false;
        break;
      }
    }
    if (needinit){
      fprintf(stderr, "Unable to locate INAQS config file, '%s', terminating!\n",
              config_names[0]);
      exit(-1);
    }
  }

  return gifs_get_forces(qm_crd, NULL, nmm, mm_crd,  mm_chg, f, f + nqm*3);
};
