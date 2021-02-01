#ifndef __GIFS_DECOHERENCE_HPP
#define __GIFS_DECOHERENCE_HPP

#include "qm_interface.hpp"
#include "electronic.hpp"

/*
  Random numbers are to be generated via the armadillo interface,
  which is guranteed to be appropriately seeded before control passes
  to Decoherence.
*/
class Decoherence {
public:
  explicit Decoherence (QMInterface *qm, double dtc): qm{qm}, dtc{dtc} {};
  virtual ~Decoherence(void) {};

  /*
    FIXME: states in Electronic are indexed starting from min_state,
    but Decoherence doesn't know what min_state is and therefore
    cannot make calls to the QMInterface.
  */
  virtual bool decohere(size_t active, Electronic &c) = 0;
  virtual void hop(Electronic &c) = 0;
  
protected:
  QMInterface *qm;
  double dtc;
  std::mt19937_64 mt64_generator;
};

#endif
