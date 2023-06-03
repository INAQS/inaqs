#ifndef GIFS_DIABATIC_SEAM_DYNAMICS_HPP
#define GIFS_DIABATIC_SEAM_DYNAMICS_HPP

#include <armadillo>
#include "configreader.hpp"
#include "bomd.hpp"


class DiabaticSeam: public BOMD{
public:
  explicit DiabaticSeam(std::shared_ptr<INAQSShared> shared, arma::mat& qm_grd, arma::mat& mm_grd): BOMD(shared, qm_grd, mm_grd) {}
  ~DiabaticSeam() {}
  double update_gradient(void) override;

protected:
  // Called within BOMD initialization
  void get_reader_data(ConfigBlockReader& reader) override; 
  ConfigBlockReader setup_reader() override;

  arma::cube g_qm, g_mm;
  arma::cube gd_qm, gd_mm;
  
  arma::mat H;  // diabatic Hamiltonian
  bool single_diabat = false;
  int which_diabat = -1;
  
  // FIXME: should probaly regularize all code on uwords etc.
  size_t upper, lower; // indicies of the adiabats to diabatize
};


#endif
