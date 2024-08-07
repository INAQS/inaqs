#include "eshake.hpp"
#include "gmx_fatal.h"

real inaqs_eshake(int natoms, real * invmass_ptr, int econq, real * qs, real invdt){
  (void) natoms;
  (void) invmass_ptr;
  (void) econq;
  (void) qs;
  (void) invdt;
  
  gmx_warning("inaqs_eshake stub called; skipping constraint step of type %d", econq);
  return 0;
}
