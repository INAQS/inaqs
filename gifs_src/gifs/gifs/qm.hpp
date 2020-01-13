#ifndef GIFS_SH_QM_CORE_H
#define GIFS_SH_QM_CORE_H

#include <vector>


class QMInterface{
public:
  QMInterface(size_t nqm, std::vector<int> &qmid);

  void get_properties(PropMap &props) {};

  template<typename T>
  void update(const T* crdqm, size_t nmm, const T* crdmm, const T* chgmm);

  inline size_t get_nqm() const noexcept { return NQM; };

protected:
  //  fixed size
  size_t NQM = 0;
  std::vector<int> atomids{};        // NQM
  std::vector<double> crd_qm{};      // NQM*3
  // flexible
  size_t NMM = 0;
  std::vector<double> crd_mm{};      // NMM*3
  std::vector<double> chg_mm{};      // NMM
  bool first_call = true;
};

#endif
