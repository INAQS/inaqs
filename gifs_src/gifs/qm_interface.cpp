#include <vector>
#include "properties.hpp"
#include "qm_interface.hpp"

QMInterface::QMInterface(const std::vector<int> &qmid, int charge, int mult):
  qm_charge(charge), qm_multiplicity(mult){
  NQM = qmid.size();
  NMM = 0;
  atomids = qmid;
  
  crd_qm.resize(3 * NQM);
}

void QMInterface::update(const std::vector<double> &crdqm, const std::vector<double> &crdmm, const std::vector<double> &chgmm){
  NMM = chgmm.size();
  crd_qm = crdqm;
  crd_mm = crdmm;
  chg_mm = chgmm;
}
