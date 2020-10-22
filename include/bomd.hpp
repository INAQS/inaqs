#ifndef GIFS_SH_BOMD_CORE_H
#define GIFS_SH_BOMD_CORE_H


#include <armadillo>
#include "qm_interface.hpp"
#include "configreader.hpp"


class BOMD
{
public:
  explicit BOMD(arma::mat& qm_grd, arma::mat& mm_grd);
  virtual ~BOMD() {delete qm;}
  
  virtual double update_gradient(void);
  virtual bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double total_energy);
  // setup from user input
  void setup(FileHandle& fh,
             const arma::uvec& atomicnumbers,
	         arma::mat& qm_crd,
	         arma::mat& mm_crd,
	         arma::vec& mm_chg
          );
  //
protected:
  virtual ConfigBlockReader setup_reader(void);
  virtual void get_reader_data(ConfigBlockReader& reader);
  void add_necessary_keys(ConfigBlockReader& reader);
  inline arma::uword NQM(void) const { return qm_grd.n_cols; }
  inline arma::uword NMM(void) const { return mm_grd.n_cols; }
  
  QMInterface* qm{nullptr};
  // fixed size
  arma::mat& qm_grd;
  // flexible size
  arma::mat& mm_grd;
  arma::vec energy{};

  double elast = 0, edrift = 0;
};

class PrintBomd:
    public BOMD 
{
public:
  explicit PrintBomd(arma::mat& qm_grd,
		             arma::mat& mm_grd) : BOMD(qm_grd, mm_grd) {}
  
  double update_gradient(void) {
      qm->crd_qm.print("qm_crd: ");
      qm->crd_mm.print("mm_crd: ");
      qm_grd.fill(0.0);
      mm_grd.fill(0.0);
      return 0.0;
  }
    //
  bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift) {
    (void) e_drift;
      velocities.print("Velocities:");
      masses.print("Masses:");
      total_gradient.print("Total Gradient");
      return true;
  };
    
private:
  virtual ConfigBlockReader setup_reader() {
    return ConfigBlockReader{"printbomd"};    
  };
  //virtual void get_reader_data(ConfigBlockReader& reader) {}

};


#endif
