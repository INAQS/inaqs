#include "properties.hpp"
#include "qm_interface.hpp"

QMInterface::QMInterface(arma::uvec& in_qmids, 
                         arma::mat& in_qm_crd, 
                         arma::mat& in_mm_crd, 
                         arma::vec& in_mm_chg, 
                         int charge, 
                         int mult,
			 int excited_states) :
  NQM {in_qmids.size()},
  qm_charge(charge), qm_multiplicity(mult),
  excited_states(excited_states),
  atomids{in_qmids}, crd_qm{in_qm_crd}, 
  NMM {in_mm_chg.size()}, crd_mm{in_mm_crd}, chg_mm{in_mm_chg}
{}


void QMInterface::update()
{
  NMM = chg_mm.size();
  qm_call_idx++;
}
