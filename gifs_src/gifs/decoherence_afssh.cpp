#include "decoherence_afssh.hpp"

/* This is steps */
bool AFSSH::decohere(size_t active, Electronic &c){
  (void) c;
  std::cout << "AFSSH::decohere() on state " << active << std::endl;
  return false;
}


/* Jain 2016 step 6*/
void AFSSH::hop(Electronic &c){
  (void) c;
  //moments.reset();
}
