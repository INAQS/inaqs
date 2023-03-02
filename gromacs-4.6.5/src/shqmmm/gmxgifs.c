#include <stdio.h>
#include "gmxgifs.h"
#include "gifs.hpp"
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>

void gifs_scale_velocities(real energy, rvec *v, rvec *f, real* invmass) {
  _Generic(
           energy,
           float: gifs_rescale_velocities_float,
           double: gifs_rescale_velocities_double)(energy, (real*) f, invmass, (real*) v);
};

void gifs_update_global_idx(int* indexQM, int* indexMM) {
  gifs_update_global_index(indexQM, indexMM);
};

real gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const rvec * qm_crd,
			size_t nmm, const rvec* mm_crd, const real* mm_chg,
			rvec* f, rvec* fshift){

  if (!gifs_interface_is_ready()){
    fprintf(stderr, "[INAQS] WARNING: attempt to call %s before init().\n", __func__);
  }

  real energy =  _Generic(
    qm_crd,
    const float (*)[3]: gifs_get_forces_float,
    const double (*)[3]: gifs_get_forces_double)((real*) qm_crd, NULL, nmm,
                                                 (real*) mm_crd,  mm_chg,
                                                 (real*) f,
                                                 (real*) (f + nqm));
  return energy;
};

bool inaqs_init(char * inaqsConfigFile, real timeStep, size_t nqm, const int * qm_atomids){
  (void) inaqsConfigFile;  // FIXME: need to actually use
  (void) timeStep;         // FIXME: need to actually use
  
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
