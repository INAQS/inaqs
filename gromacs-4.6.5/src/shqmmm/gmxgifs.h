#ifndef GMX_GIFS_INTERFACE_H
#define GMX_GIFS_INTERFACE_H
#include <stddef.h>
#include <stdbool.h>

void gifs_update_global_idx(int* indexQM, int* indexMM);
void gifs_scale_velocities(float energy, float *v, float *f, float* invmass);
float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
			size_t nmm, const float* mm_crd, const float* mm_chg,
                       float* f, float* fshift);

bool inaqs_init(char * configfile, float timeStep, size_t nqm, const int * qm_atomids);

/*

FIXME: Compilation errors to be dealt with.

/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c: In function ‘call_QMroutine’:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:170:67: warning: passing argument 3 of ‘gifs_do_qm_forces’ from incompatible pointer type [-Wincompatible-pointer-types]
     QMener = gifs_do_qm_forces(qm->nrQMatoms, qm->atomicnumberQM, qm->xQM,
                                                                   ^~
In file included from /data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:45:0:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/../shqmmm/gmxgifs.h:6:7: note: expected ‘const float *’ but argument is of type ‘real (*)[3] {aka float (*)[3]}’
 float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
       ^~~~~~~~~~~~~~~~~
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:171:26: warning: passing argument 5 of ‘gifs_do_qm_forces’ from incompatible pointer type [-Wincompatible-pointer-types]
           mm->nrMMatoms, mm->xMM, mm->MMcharges,
                          ^~
In file included from /data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:45:0:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/../shqmmm/gmxgifs.h:6:7: note: expected ‘const float *’ but argument is of type ‘real (*)[3] {aka float (*)[3]}’
 float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
       ^~~~~~~~~~~~~~~~~
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:172:11: warning: passing argument 7 of ‘gifs_do_qm_forces’ from incompatible pointer type [-Wincompatible-pointer-types]
           f, fshift);
           ^
In file included from /data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:45:0:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/../shqmmm/gmxgifs.h:6:7: note: expected ‘float *’ but argument is of type ‘real (*)[3] {aka float (*)[3]}’
 float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
       ^~~~~~~~~~~~~~~~~
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:172:14: warning: passing argument 8 of ‘gifs_do_qm_forces’ from incompatible pointer type [-Wincompatible-pointer-types]
           f, fshift);
              ^~~~~~
In file included from /data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/qmmm.c:45:0:
/data/home/vale/projects/GQSH/gifs/gromacs-4.6.5/src/mdlib/../shqmmm/gmxgifs.h:6:7: note: expected ‘float *’ but argument is of type ‘real (*)[3] {aka float (*)[3]}’
 float gifs_do_qm_forces(size_t nqm, const int* qm_atomids, const float* qm_crd,
       ^~~~~~~~~~~~~~~~~

*/

#endif
