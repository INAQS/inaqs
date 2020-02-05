#ifndef GIFS_SH_QM_INTERFACE_HPP
#define GIFS_SH_QM_INTERFACE_HPP

#include "properties.hpp"
#include <vector>

/*
  Abstract base class for all implementations of QMInterface 
*/

class QMInterface{
public:
  QMInterface(std::vector<int>& qmids, std::vector<double>& qm_crd, std::vector<double>& mm_crd, std::vector<double>& mm_chg, int charge, int mult);
  void update();
  inline size_t nqm() const noexcept { return NQM; };
  inline int call_idx() const noexcept { return qm_call_idx; };
  
  virtual void get_properties(PropMap &props) =0;
  virtual ~QMInterface(){};

protected:
  //Properties
  size_t NQM;             // const, actually NQM+NLink
  int qm_charge, qm_multiplicity;
  //  fixed size
  std::vector<int>& atomids;        // NQM
  std::vector<double>& crd_qm;      // NQM*3
  // flexible
  size_t NMM;
  std::vector<double>& crd_mm;      // NMM*3
  std::vector<double>& chg_mm;      // NMM

private:
  int qm_call_idx = 0; // tracks each call to update;
};

#endif
