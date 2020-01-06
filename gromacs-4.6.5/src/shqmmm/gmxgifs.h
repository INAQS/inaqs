#ifndef GMX_GIFS_INTERFACE_H
#define GMX_GIFS_INTERFACE_H
#include <stddef.h>

void gifs_scale_velocities(float* v, float* f, float* invmass);
float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
			size_t nmm, const float* mm_crd, const float* mm_chg,
                       float* f, float* fshift);

/*

FIXME: Compilation errors to be dealt with.

/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/kernel/md.c: In function ‘do_md’:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/kernel/md.c:1246:35: warning: passing argument 1 of ‘gifs_scale_velocities’ from incompatible pointer type [-Wincompatible-pointer-types]
             gifs_scale_velocities(state->v, f, mdatoms->invmass);
                                   ^~~~~
In file included from /data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/kernel/md.c:108:0:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/kernel/../shqmmm/gmxgifs.h:5:6: note: expected ‘float *’ but argument is of type ‘real (*)[3] {aka float (*)[3]}’
 void gifs_scale_velocities(float* v, float* f, float* invmass);
      ^~~~~~~~~~~~~~~~~~~~~
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/kernel/md.c:1246:45: warning: passing argument 2 of ‘gifs_scale_velocities’ from incompatible pointer type [-Wincompatible-pointer-types]
             gifs_scale_velocities(state->v, f, mdatoms->invmass);
                                             ^
In file included from /data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/kernel/md.c:108:0:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/kernel/../shqmmm/gmxgifs.h:5:6: note: expected ‘float *’ but argument is of type ‘real (*)[3] {aka float (*)[3]}’
 void gifs_scale_velocities(float* v, float* f, float* invmass);

*/

#endif
