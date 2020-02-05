#include <vector>
#include "properties.hpp"
#include "qm_interface.hpp"

QMInterface::QMInterface(std::vector<int>& in_qmids, 
                         std::vector<double>& in_qm_crd, 
                         std::vector<double>& in_mm_crd, 
                         std::vector<double>& in_mm_chg, 
                         int charge, 
                         int mult) :
    atomids{in_qmids}, crd_qm{in_qm_crd}, 
    crd_mm{in_mm_crd}, chg_mm{in_mm_chg}, 
    qm_charge(charge), qm_multiplicity(mult)
{
  NQM = atomids.size();
  NMM = 0;
  
  crd_qm.resize(3 * NQM);
};


void QMInterface::update()
{
  NMM = chg_mm.size();
  qm_call_idx++;
}
