#ifndef __QM_HPP
#define __QM_HPP

#include "properties.hpp"
#include <vector>

class QMInterface{
public:
  QMInterface(size_t nqm, std::vector<int> &qmid);
  void get_properties(PropMap &props);

  void update(std::vector<double> &crdqm,
	      std::vector<double> &crdmm,
	      std::vector<double> &chgmm);  

  void update(const float* crdqm, size_t nmm, const float* crdmm, const float* chgmm);
  inline size_t get_nqm() const noexcept { return NQM; };

protected:
  void get_gradient_energies(std::vector<double> &g_qm, std::vector<double> &g_mm, std::vector<double> &e);
  void write_gradient_job(std::string fname);
  void exec_qchem(std::string qcprog, std::string ifname, std::string savdir);
  void parse_qm_gradient(std::string savdir, std::vector<double> &g_qm, std::vector<double> &e);
  std::string get_qcprog(void);
  size_t NQM;             // const, actually NQM+NLink
  //  fixed size
  std::vector<int> atomids;        // NQM
  std::vector<double> crd_qm;      // NQM*3
  // flexible
  size_t NMM;
  std::vector<double> crd_mm;      // NMM*3
  std::vector<double> chg_mm;      // NMM
  bool first_call = true;
};

#endif //__QM_HPP
