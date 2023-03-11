#ifndef GIFS_DIABATIC_SEAM_DYNAMICS_HPP
#define GIFS_DIABATIC_SEAM_DYNAMICS_HPP

#include <armadillo>
#include "configreader.hpp"
#include "bomd.hpp"


class DiabaticSeam: public BOMD{
public:
  explicit DiabaticSeam(double classicalTimeStep, arma::mat& qm_grd, arma::mat& mm_grd): BOMD(classicalTimeStep, qm_grd, mm_grd) {}
  ~DiabaticSeam() {}
  double update_gradient(void) override;

protected:
  // Called within BOMD initialization
  void get_reader_data(ConfigBlockReader& reader) override; 
  ConfigBlockReader setup_reader() override;
  
  double build_diabatic_forces_projected(void);
  double build_diabatic_forces_restrained(void);

  arma::cube g_qm, g_mm;
  arma::cube gd_qm, gd_mm;
  
  arma::vec diabatic_energy;  
  arma::mat diabatic_rot_mat;

  double alpha = 0;

  bool single_diabat = false;
  
  // FIXME: should probaly regularize all code on uwords etc.
  size_t upper, lower; // indicies of the adiabats to diabatize
};


#endif
