#ifndef GMX_GIFS_INTERFACE_H
#define GMX_GIFS_INTERFACE_H

void gifs_say_hallo(int fp);
void gifs_scale_velocities(float* v, float* f, float* invmass);
void gifs_do_qm_forces(const int nqm, const int* qm_atomids, const float* qm_crd, 
                       const int nmm, const float* mm_crd, const float* mm_chg,
                       float* f, float* fshift, 
                       const float BOHR2NM, const float HARTREE_BOHR2MD);
#endif
