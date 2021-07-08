#ifndef GIFS_SH_QM_INTERFACE_HPP
#define GIFS_SH_QM_INTERFACE_HPP

#include "properties.hpp"
#include <armadillo>

/*
  Abstract base class for all implementations of QMInterface 
*/

class QMInterface{
  friend class PrintBomd;
public:
  QMInterface(const arma::uvec& in_qmids, 
	      arma::mat& in_qm_crd, 
	      arma::mat& in_mm_crd, 
	      arma::vec& in_mm_chg, 
	      const int charge, 
	      const int mult,
	      const int excited_states,
              const int min_state);
  virtual void update(); // if overriding, be sure to call the parent too.
  inline size_t nqm() const noexcept { return NQM; };
  inline int call_idx() const noexcept { return qm_call_idx; };
  
  virtual void get_properties(PropMap &props) =0;
  virtual ~QMInterface(){};

protected:
  //Properties
  const size_t NQM;             // const, actually NQM+NLink
  const int qm_charge;
  const int qm_multiplicity;
  size_t excited_states;        // non-const since a child might want to override
  const size_t min_state;
  
  //  fixed size
  const arma::uvec& atomids;    // NQM
  arma::mat& crd_qm;      // NQM*3
  // flexible
  size_t NMM;
  arma::mat& crd_mm;      // NMM*3
  arma::vec& chg_mm;      // NMM

private:
  int qm_call_idx = 0; // tracks each call to update;
};

#endif
