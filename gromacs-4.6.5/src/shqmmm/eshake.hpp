#ifndef __INAQS_ESHAKE_HPP
#define __INAQS_ESHAKE_HPP

#include "typedefs.h"
#include "nrnb.h"
#include "constr.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  real inaqs_eshake(int natoms, real * invmass_ptr, int econq, real * qs, real invdt);  
  /* Shake steps for electronic constraints*/

#ifdef __cplusplus
}
#endif
#endif
