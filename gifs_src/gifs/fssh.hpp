#ifndef __GIFS_FSSH_HPP
#define __GIFS_FSSH_HPP

#include "bomd.hpp"
#include <armadillo>

class FSSH: public BOMD{
public:
  explicit FSSH(int nqm, const int * qmid); // need to parse our config
  virtual ~FSSH() {};
  
protected:
  void main(void);
  double gen_rand(void);

  arma::cx_vec c;
  arma::mat U;
  arma::mat T, V;

  double dtc;
  double dtq;
  
  size_t min_state;
  size_t excited_states;
  size_t active_state;
  size_t target_state;
  bool hopping = false;

  std::mt19937_64 mt64_generator;
};

#endif
