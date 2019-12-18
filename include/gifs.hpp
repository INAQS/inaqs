#ifndef GIFS_MAIN_H
#define GIFS_MAIN_H

/*Constants from gromacs*/
#define AVOGADRO         (6.0221367e23)
#define HARTREE2KJ       (4.3597482e-21)
#define BOHR2NM          (0.0529177249)
#define HARTREE_BOHR2MD  (HARTREE2KJ*AVOGADRO/BOHR2NM)

#ifdef __cplusplus
extern "C" {
#endif
 
  void create_qm_interface(size_t nqm, const int* qm_atomids);
  float gifs_get_forces(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift);

#ifdef __cplusplus 
}
#endif

#endif
