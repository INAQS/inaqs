#ifndef GIFS_QM_TULLY_HPP
#define GIFS_QM_TULLY_HPP

#include "properties.hpp"
#include "qm_interface.hpp"
#include "configreader.hpp"
#include <armadillo>
#include <memory>

class QMModel{
public:
  QMModel(void);

  virtual void update(const arma::vec &x) =0;
  
  const arma::vec &              energies(void)  const {return energy;}
  const arma::cube &             gradients(void) const {return gradient;}
  const arma::field<arma::mat> & nacs(void)      const {return nac;}
  const arma::mat &              overlaps(void)  const {return overlap;}
  
protected:
  arma::vec energy;
  arma::cube gradient;
  arma::field<arma::mat> nac;
  arma::mat overlap;
};


class Tully1: public QMModel{
public:
  Tully1(void);
  void update(const arma::vec &x) override {(void) x;};
};


class QMTully: public QMInterface{
public:
  QMTully(FileHandle& fh, 
           const arma::uvec& in_qmids, 
	   arma::mat& in_qm_crd, 
	   arma::mat& in_mm_crd, 
	   arma::vec& in_mm_chg, 
	   const int charge, 
	   const int mult,
	   const int excited_states,
           const int min_state);
  
  void get_properties(PropMap &props) override;
  void update(void) override;

private:
  ConfigBlockReader tully_reader();
  std::unique_ptr<QMModel> model;
};

#endif
