#ifndef GIFS_SH_QM_QCHEM_HPP
#define GIFS_SH_QM_QCHEM_HPP

#include "properties.hpp"
#include "qm_interface.hpp"
#include <vector>

class QM_QChem: public QMInterface{
public:
  void get_properties(PropMap &props);
  QM_QChem(const std::vector<int> &qmid, int charge, int mult);
  
private:
  inline void ang2bohr(std::vector<double> &v);
  void get_gradient_energies(std::vector<double> &g_qm, std::vector<double> &g_mm, std::vector<double> &e);
  void write_gradient_job(void);
  void write_molecule_section(std::ostream &ifile);
  void exec_qchem(void);
  void parse_qm_gradient(std::vector<double> &g_qm, std::vector<double> &e);
  void parse_mm_gradient(std::vector<double> &g_mm);
  std::string get_qcprog(void);
  std::string get_qcscratch(void);
  int readQFMan(int filenum, std::vector<double> &v);

  const std::string qc_scratch_directory;
  const std::string qc_executable;
  const std::string qc_input_file = "GQSH.in";
  bool first_call = true;
};

#endif
