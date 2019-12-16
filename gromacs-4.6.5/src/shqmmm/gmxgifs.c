#include <stdio.h>
#include "gmxgifs.h"

void
gifs_say_hallo(int fp)
{
    fprintf(fp, "HALLO");
}


void gifs_scale_velocities(float* v, float* f, float* invmass) { }

float gifs_do_qm_forces(const int nqm, const int* qm_atomids, const float* qm_crd, 
                       const int nmm, const float* mm_crd, const float* mm_chg,
                       float* f, float* fshift)
{
    return 0.0;
}

