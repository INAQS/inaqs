#ifndef __GIFS_DECOHERENCE_HPP
#define __GIFS_DECOHERENCE_HPP

#include "qm_interface.hpp"
#include "electronic.hpp"

/*
  N.B.: Random numbers are to be generated via the armadillo
  interface, which is guaranteed to be appropriately seeded before
  control passes to Decoherence via decohere() or hopped().
*/
class Decoherence {
public:
  explicit Decoherence (QMInterface ** const qm, const double dtc,
                        const size_t min_state,
                        const size_t hopping_states,
                        const size_t nqm, const size_t nmm):
    qm{qm}, dtc{dtc}, min_state{min_state}, nstates{hopping_states}, nqm{nqm}, nmm{nmm} {};
  virtual ~Decoherence(void) {};

  virtual bool decohere(Electronic &c, const arma::mat U, const size_t active_state, const arma::vec v, const arma::vec m) = 0;
  virtual void hopped(Electronic &c, size_t active_state) = 0;
  
protected:
  QMInterface ** const qm;
  const double dtc;
  const size_t min_state;
  const size_t nstates;
  const size_t nqm;
  const size_t nmm;
};

#endif
