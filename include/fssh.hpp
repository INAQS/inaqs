#ifndef __GIFS_FSSH_HPP
#define __GIFS_FSSH_HPP

#include <armadillo>
//
#include "bomd.hpp"
#include "configreader.hpp"
#include "electronic.hpp"

class FSSH: public BOMD{
public:
  //explicit FSSH(int nqm, const int * qmid, size_t min_state, size_t excited_states, size_t active_state, double dtc);
  explicit FSSH(
        FileHandle& fh,
        arma::uvec& atomicnumbers,  // need to parse our config
		arma::mat& qm_crd, 
		arma::mat& mm_crd, 
		arma::vec& mm_chg, 
		arma::mat& qm_grd,
		arma::mat& mm_grd);
  virtual ~FSSH() {};

  virtual double update_gradient(void);
  virtual bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift);
  
protected:
  virtual void get_reader_data(ConfigBlockReader& reader); 
  virtual ConfigBlockReader setup_reader();
  void electonic_evolution(void);
  void update_md_global_gradient(arma::mat &total_gradient);
  void hop_and_scale(arma::mat &velocities, arma::vec &mass);
  double gen_rand(void);
  arma::uword sample_discrete(const arma::vec &p);

  arma::mat U;
  arma::mat T, V;

  double dtc;
  double delta_e_tol;
  double dtq;

  size_t min_state;
  size_t excited_states;
  
  size_t active_state;
  size_t target_state;
  bool hopping = false;

  Electronic c {};

  arma::mat nac;
  
  std::mt19937_64 mt64_generator;
};

#endif
