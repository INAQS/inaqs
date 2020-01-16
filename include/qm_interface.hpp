#ifndef GIFS_SH_QM_INTERFACE_HPP
#define GIFS_SH_QM_INTERFACE_HPP

#include "properties.hpp"
#include <vector>

/*
  Virtual base class for all implementations of QMInterface 
*/

class QMInterface{
public:
  QMInterface(const std::vector<int> &qmid, int charge, int mult);
  void update(const std::vector<double> &crdqm, const std::vector<double> &crdmm, const std::vector<double> &chgmm);
  inline size_t nqm() const noexcept { return NQM; };
  
  virtual void get_properties(PropMap &props) =0;
  virtual ~QMInterface(){};

protected:
  //Properties
  size_t NQM;             // const, actually NQM+NLink
  int qm_charge, qm_multiplicity;
  //  fixed size
  std::vector<int> atomids;        // NQM
  std::vector<double> crd_qm;      // NQM*3
  // flexible
  size_t NMM;
  std::vector<double> crd_mm;      // NMM*3
  std::vector<double> chg_mm;      // NMM
};

#endif
