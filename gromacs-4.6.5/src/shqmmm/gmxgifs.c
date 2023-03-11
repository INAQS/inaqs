#include <stdio.h>
#include "gmxgifs.h"
#include "gifs.hpp"
#include "typedefs.h"
#include "gmx_fatal.h"
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

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
    gmx_fatal(FARGS, "Attempt to use INAQS before init()");
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

/*
  FIXME: support restarts
  FIXME: consider wraping all the MD parameters (dt, restart...) in a struct
*/
bool inaqs_init(const char * inaqsConfigFile, real classicalTimeStep, size_t nqm, const int * qm_atomids){
  if (inaqsConfigFile){
    create_qm_interface(inaqsConfigFile, classicalTimeStep, nqm, qm_atomids);
    return true;
  }
  
  const char * config_names[] = {
    "inaqs_config.dat", // default name is first
    "inaqs_config.ini", // legacy inputs later
    "gifs_config.ini",
    NULL                // NULL-terminate
  };

  gmx_warning("No config file passed to mdrun via -inaqs; attempting to find the default: %s", config_names[0]);
  
  for (const char ** fname = config_names; *fname; fname++){
    if (0 == access(*fname, F_OK)){
      if (fname != config_names){
        gmx_warning("DEPRECATION NOTICE, you are not using the standard config file"
                    ", '%s'. Please update your inputs.\n", config_names[0]);
      }
      create_qm_interface(*fname, classicalTimeStep, nqm, qm_atomids);
      return true;
    }
  }

  gmx_fatal(FARGS, "Unable to locate INAQS config file, %s", config_names[0]);
  exit(-1);
}
