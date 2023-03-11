#ifndef INAQS_EHRENFEST_HPP
#define INAQS_EHRENFEST_HPP

#include <armadillo>
#include "bomd.hpp"
#include "electronic.hpp"
#include "configreader.hpp"

class Ehrenfest: public BOMD {
public:
  explicit Ehrenfest(double classicalTimeStep, arma::mat& qm_grd, arma::mat& mm_grd): BOMD{classicalTimeStep, qm_grd, mm_grd} {}
  ~Ehrenfest() { }
  double update_gradient(void) override;

protected:
  // Called within BOMD initialization
  void get_reader_data(ConfigBlockReader& reader) override; 
  ConfigBlockReader setup_reader() override;
  
  arma::mat U, T, V;
  arma::vec phases;
  arma::cube qm_grds, mm_grds;

  double dtq;
  
  size_t min_state;
  int nstates;

  std::string amplitude_file;
  
  // size_t active_state; inherited from BOMD
  Electronic c {};
};

#endif
