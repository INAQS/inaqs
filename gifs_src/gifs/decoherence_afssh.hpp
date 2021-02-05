#ifndef __GIFS_DECOHERENCE_AFSSH_HPP
#define __GIFS_DECOHERENCE_AFSSH_HPP
#include "decoherence.hpp"

/*
  Implements decoherence from the efficient "Augmented" FSSH scheme
  outlined in Jain et al. JCTC 2016.
*/

class AFSSH: public Decoherence{
public:
  explicit AFSSH(QMInterface ** const qm, const double dtc,
                 const size_t min_state,
                 const size_t shstates,
                 const size_t nqm, const size_t nmm);
  virtual ~AFSSH() {};

  bool decohere(Electronic &c, const arma::mat U, const size_t active_state, const arma::vec v, const arma::vec m) override;
  void hopped(Electronic &c, size_t active_state) override;

private:
  void evolve_moments(const Electronic &c, const arma::mat U, const size_t a, const arma::vec m);
  void build_invtau(arma::vec &invtau_d, arma::vec &invtau_r, const arma::mat U, const size_t a, const arma::vec v);
  
  void collapse(Electronic &c, size_t active_state, size_t collapse_state);
  void reset_moments(void);

  // current energies
  arma::vec V;
  
  // adiabatic moments
  arma::field<arma::vec> dR, dP, dF;
  // diabatic moments
  arma::field<arma::vec> dR_, dP_, dF_;
};

#endif
