#ifndef GIFS_SH_BOMD_RESCALE_CORE_H
#define GIFS_SH_BOMD_RESCALE_CORE_H

#include <armadillo>
#include "bomd.hpp"


class RescaleBomd: 
    public BOMD 
{
public:
  explicit RescaleBomd(arma::mat& qm_grd, arma::mat& mm_grd): BOMD(qm_grd, mm_grd) {}
      
  ~RescaleBomd() {delete qm;}
  
  virtual double update_gradient(void) override;
  virtual bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift) override;
    //
protected:
  virtual ConfigBlockReader setup_reader(void) override;
  virtual void get_reader_data(ConfigBlockReader& reader) override;
private:
  double additional_energy{100.0};
  double dE{0.01};
  int step{10};
  double total_energy{0.0};
};



#endif
