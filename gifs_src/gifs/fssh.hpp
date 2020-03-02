#ifndef __GIFS_FSSH_HPP
#define __GIFS_FSSH_HPP

#include "bomd.hpp"
#include "electronic.hpp"
#include <armadillo>

class FSSH: public BOMD{
public:
  explicit FSSH(int nqm, const int * qmid, size_t min_state, size_t excited_states, size_t active_state, double dtc);
  virtual ~FSSH() {};
  
protected:
  void electonic_evolution(void);
  void hop_and_scale(arma::vec vel, arma::vec inv_mass);
  double gen_rand(void);

  arma::mat U;
  arma::mat T, V;

  const double dtc;
  double dtq;

  const size_t min_state;
  const size_t excited_states;
  size_t active_state;
  size_t target_state;
  bool hopping = false;

  Electronic c;

  arma::mat nac;
  
  std::mt19937_64 mt64_generator;
};

#endif
