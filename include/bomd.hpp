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
  inline arma::uword NQM(void) const { return qm_grd.n_cols; }
  inline arma::uword NMM(void) const { return mm_grd.n_cols; }
  inline int call_idx() const noexcept { return md_call_idx; };

  template <typename T>
  void saveh5(const T &tensor, std::string kind){
    tensor.save(
                arma::hdf5_name("gifs.hdf5", "/" + kind + "/" +
                                std::to_string(call_idx()),
                                arma::hdf5_opts::replace)
                );
  };

  /* Config for keys common to all dynamics classes */
  void add_common_keys(ConfigBlockReader& reader);
  
  /* Child classes will override these methods for own setup */
  virtual ConfigBlockReader setup_reader(void);
  virtual void get_reader_data(ConfigBlockReader& reader);
  
  QMInterface* qm{nullptr};
  // fixed size
  arma::mat& qm_grd;
  // flexible size
  arma::mat& mm_grd;
  arma::vec energy{};

  // On which surface are we running?
  size_t active_state;

  // To track energy drift from GMX
  double elast = 0, edrift = 0;
  
private:
  //FIXME: need to unify _call_idx() tracking at a higher level (gifs?)
  int md_call_idx = 1; // tracks each call to rescale_velocities();
};

class PrintBomd:
    public BOMD 
{
public:
  virtual ~PrintBomd() {};
  explicit PrintBomd(arma::mat& qm_grd,
                     arma::mat& mm_grd) : BOMD(qm_grd, mm_grd) {}
  
  double update_gradient(void) {
    double e = BOMD::update_gradient();
    qm->crd_qm.t().print("qm_crd: ");
    qm->crd_mm.t().print("mm_crd: ");
    qm_grd.t().print("qm_grd");
    mm_grd.t().print("mm_grd");
    return e;
  }

  bool rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double total_energy) {
    BOMD::rescale_velocities(velocities, masses, total_gradient, total_energy);
    velocities.t().print("Velocities:");
    masses.t().print("Masses:");
    total_gradient.t().print("Total Gradient");
    return true;
  }
};


#endif
