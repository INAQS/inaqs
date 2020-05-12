#ifndef GIFS_CONSTANTS_H
#define GIFS_CONSTANTS_H

#define AU2SI_LEN  5.29177210903e-11;   // Bohr: meters     
#define AU2SI_MASS 9.1093837015e-31;    // kilogram: au     
#define AU2SI_TIME 2.418884326509e-17;  // seconds:  au

// From Gromacs
#define AVOGADRO         (6.0221367e23)
#define HARTREE2KJ       (4.3597482e-21)
#define BOHR2NM          (0.0529177249)
#define HARTREE_BOHR2MD  (HARTREE2KJ*AVOGADRO/BOHR2NM)

#endif
