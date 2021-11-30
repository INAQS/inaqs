#ifndef GIFS_BOMD_ELECTRONIC_HPP
#define GIFS_BOMD_ELECTRONIC_HPP

#include <armadillo>
#include "bomd.hpp"
#include "electronic.hpp"
#include "configreader.hpp"

class ElectronicBomd: public BOMD {
public:
  explicit ElectronicBomd(arma::mat& qm_grd, arma::mat& mm_grd): BOMD{qm_grd, mm_grd} {}
  ~ElectronicBomd() { }
  double update_gradient(void) override;

protected:
  // Called within BOMD initialization
  void get_reader_data(ConfigBlockReader& reader) override; 
  ConfigBlockReader setup_reader() override;
  
  arma::mat U, T, V;

  double dtc;
  double dtq;
  
  size_t min_state;
  int nstates;

  std::string amplitude_file;
  
  // size_t active_state; inherited from BOMD
  Electronic c {};
};

#endif
