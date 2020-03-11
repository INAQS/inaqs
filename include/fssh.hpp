#ifndef __GIFS_FSSH_HPP
#define __GIFS_FSSH_HPP

#include "bomd.hpp"
#include "electronic.hpp"
#include <armadillo>

class FSSH: public BOMD{
public:
  //explicit FSSH(int nqm, const int * qmid, size_t min_state, size_t excited_states, size_t active_state, double dtc);
  explicit FSSH(arma::uvec& atomicnumbers,  // need to parse our config
		arma::mat& qm_crd, 
		arma::mat& mm_crd, 
		arma::vec& mm_chg, 
		arma::mat& qm_grd,
		arma::mat& mm_grd, size_t min_state, size_t excited_states, size_t active_state, double dtc,  double delta_e_tol);
  virtual ~FSSH() {};

  virtual double update_gradient(void);
  virtual void rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift);
  
protected:
  void electonic_evolution(void);
  void update_md_global_gradient(arma::mat &total_gradient);
  void hop_and_scale(arma::mat &velocities, arma::vec &mass);
  double gen_rand(void);
  arma::uword sample_discrete(const arma::vec &p);

  arma::mat U;
  arma::mat T, V;

  const double dtc;
  const double delta_e_tol;
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
