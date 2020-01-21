#ifndef GIFS_SH_QM_QCHEM_HPP
#define GIFS_SH_QM_QCHEM_HPP

#include "properties.hpp"
#include "qm_interface.hpp"
#include <vector>
#include <map>

class QM_QChem: public QMInterface{
public:
  QM_QChem(const std::vector<int> &qmid, int charge, int mult);
  void get_properties(PropMap &props);
  
private:
  inline void ang2bohr(std::vector<double> &v);
  void get_gradient_energies(std::vector<double> &g_qm, std::vector<double> &g_mm, std::vector<double> &e);
  void write_gradient_job(void);
  void write_molecule_section(std::ostream &ifile);
  void write_rem_section(std::ostream &os, std::map<std::string, std::string> options);
  std::ofstream get_input_handle(void);
  void exec_qchem(void);
  void parse_qm_gradient(std::vector<double> &g_qm, std::vector<double> &e);
  void parse_mm_gradient(std::vector<double> &g_mm);
  std::string get_qcprog(void);
  std::string get_qcscratch(void);
  int readQFMan(int filenum, std::vector<double> &v);

  const std::string qc_scratch_directory;
  const std::string qc_executable;
  const std::string qc_input_file = "GQSH.in";
  const std::string qc_log_file = "GQSH.out";
  const std::string exchange_method;
  const std::string basis_set;
  bool first_call = true;
};

#endif
