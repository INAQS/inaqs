#ifndef __GIFS_FSSH_HPP
#define __GIFS_FSSH_HPP

#include <armadillo>
#include <memory>
#include "decoherence.hpp"
#include "bomd.hpp"
#include "configreader.hpp"
#include "electronic.hpp"

enum class VelocityReversal{
  none = false,
  derivative_coupling,
  derivative_coupling_velocity,
  target_gradient,
  diabatic_difference,
};

std::ostream& operator<<( std::ostream& os, const VelocityReversal& v);
VelocityReversal from_string(std::string str);

class FSSH: public BOMD{
public:
  explicit FSSH(double classicalTimeStep, arma::mat& qm_grd, arma::mat& mm_grd): BOMD{classicalTimeStep, qm_grd, mm_grd} {}
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
  void check_overlap(const arma::mat& U);

  VelocityReversal velocity_reversal = VelocityReversal::derivative_coupling;
  double trivial_crossing_threshold;

  arma::mat U, T, V;
  arma::vec phases;
  arma::mat nac;

  double delta_e_tol;
  bool rescale_initial_velocities = false;
  double dtq;
  
  size_t min_state;
  int shstates;

  std::string amplitude_file;
  
  // size_t active_state; inherited from BOMD
  size_t target_state;
  bool hopping = false;

  Electronic c {};
  Decoherence *decoherence {nullptr};
};

#endif
