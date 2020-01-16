#ifndef GIFS_SH_QM_CORE_H
#define GIFS_SH_QM_CORE_H

#include "properties.hpp"
#include <vector>

class QMInterface{
public:
  QMInterface(size_t nqm, const int * qmid);
  void get_properties(PropMap &props);

  template<typename T>
  void update(const T* crdqm, size_t nmm, const T* crdmm, const T* chgmm);

  inline size_t get_nqm() const noexcept { return NQM; };

protected:
  inline void ang2bohr(std::vector<double> &v);
  void get_gradient_energies(std::vector<double> &g_qm, std::vector<double> &g_mm, std::vector<double> &e);
  void write_gradient_job(void);
  void exec_qchem(void);
  void parse_qm_gradient(std::vector<double> &g_qm, std::vector<double> &e);
  void parse_mm_gradient(std::vector<double> &g_mm);
  std::string get_qcprog(void);
  std::string get_qcscratch(void);
  int readQFMan(int filenum, std::vector<double> &v);

  //Properties
  size_t NQM;             // const, actually NQM+NLink
  int qm_charge, qm_multiplicity;
  const std::string qc_scratch_directory;
  const std::string qc_executable;
  const std::string qc_input_file = "GQSH.in";
  //  fixed size
  std::vector<int> atomids;        // NQM
  std::vector<double> crd_qm;      // NQM*3
  // flexible
  size_t NMM;
  std::vector<double> crd_mm;      // NMM*3
  std::vector<double> chg_mm;      // NMM
  bool first_call = true;
};

#endif
