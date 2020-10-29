#ifndef __GIFS_FSSH_HPP
#define __GIFS_FSSH_HPP

#include <armadillo>
#include "bomd.hpp"
#include "configreader.hpp"
#include "electronic.hpp"

class FSSH: public BOMD{
public:
  explicit FSSH(arma::mat& qm_grd,
		        arma::mat& mm_grd);
  virtual ~FSSH() {};

  virtual double update_gradient(void);
  virtual bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift);
  
protected:
  virtual void get_reader_data(ConfigBlockReader& reader); 
  virtual ConfigBlockReader setup_reader();
  void electonic_evolution(void);
  void backpropagate_gradient_velocities(arma::mat &total_gradient, arma::mat &velocities, arma::vec &masses);
  void hop_and_scale(arma::mat &velocities, arma::vec &mass);
  arma::uword sample_discrete(const arma::vec &p);

  arma::mat U;
  arma::mat T, V;

  double dtc;
  double delta_e_tol;
  double dtq;
  
  size_t min_state, excited_states;
  int shstates;

  std::string amplitude_file;
  
  // size_t active_state; inherited from BOMD
  size_t target_state;
  bool hopping = false;

  Electronic c {};

  arma::mat nac;
  
  std::mt19937_64 mt64_generator;
};

#endif
