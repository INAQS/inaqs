#ifndef __GIFS_FSSH_HPP
#define __GIFS_FSSH_HPP

#include <armadillo>
#include <memory>
#include "decoherence.hpp"
#include "bomd.hpp"
#include "configreader.hpp"
#include "electronic.hpp"


class FSSH: public BOMD{
public:
  explicit FSSH(arma::mat& qm_grd, arma::mat& mm_grd): BOMD{qm_grd, mm_grd} {}
  virtual ~FSSH() { delete decoherence; }
  virtual double update_gradient(void) override;

  //FIXME: make velocity_rescale return an energy difference
  virtual bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift) override;
  
protected:
  // Called within BOMD initialization
  virtual void get_reader_data(ConfigBlockReader& reader) override; 
  virtual ConfigBlockReader setup_reader() override;
  
  void electronic_evolution(void);
  double hop_and_scale(arma::mat &total_gradient, arma::mat &velocities, const arma::vec &m);

  arma::mat U, T, V;
  arma::mat nac;

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
  Decoherence *decoherence {nullptr};
};

#endif
