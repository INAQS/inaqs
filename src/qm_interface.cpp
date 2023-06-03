#include "properties.hpp"
#include "qm_interface.hpp"

QMInterface::QMInterface(std::shared_ptr<INAQSShared> shared,
                         const arma::uvec& in_qmids, 
                         arma::mat& in_qm_crd, 
                         arma::mat& in_mm_crd, 
                         arma::vec& in_mm_chg, 
                         const int charge, 
                         const int mult,
                         const int excited_states,
                         const int min_state) :
  shared{shared},
  NQM {in_qmids.size()},
  qm_charge(charge), qm_multiplicity(mult),
  excited_states(excited_states),
  min_state(min_state),
  atomids{in_qmids}, crd_qm{in_qm_crd}, 
  NMM {in_mm_chg.size()}, crd_mm{in_mm_crd}, chg_mm{in_mm_chg}
{
  if (excited_states < 0){
    std::cerr << "excited_states=" << excited_states << std::endl;
    throw std::logic_error("Cannot have negatively indexed state!");
  }

  if (min_state < 0){
    std::cerr << "min_state=" << min_state << std::endl;
    throw std::logic_error("Cannot have negatively indexed state!");
  }
}


void QMInterface::update()
{
  NMM = chg_mm.size();
  //qm_call_idx++;
}
