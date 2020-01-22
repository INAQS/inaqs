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
  void get_gradient_energies(std::vector<double> &g_qm, std::vector<double> &g_mm, std::vector<double> &e);
  void get_excited_gradient( std::vector<double> &g_qm, std::vector<double> &g_mm, std::vector<double> &e, size_t surface);
  void get_nac_vector(std::vector<double> &nac, size_t A, size_t B);

  void exec_qchem(void);
  void write_molecule_section(std::ostream &ifile);
  void write_rem_section(std::ostream &os, std::map<std::string, std::string> options);
  std::ofstream get_input_handle(void);

  void parse_qm_gradient(std::vector<double> &g_qm);
  void parse_energies(std::vector<double> &e);
  void parse_mm_gradient(std::vector<double> &g_mm);
  void parse_nac_vector(std::vector<double> &nac);
  size_t readQFMan(int filenum, std::vector<double> &v);
  size_t readQFMan(int filenum, std::vector<double> &v, size_t N, size_t offset);

  const std::string get_qcprog(void);
  const std::string get_qcscratch(void);

  const std::string qc_scratch_directory;
  const std::string qc_executable;
  const std::string qc_input_file = "GQSH.in";
  const std::string qc_log_file = "GQSH.out";
  const std::string exchange_method;
  const std::string basis_set;
  const size_t excited_states;
  bool first_call = true;
};

/* Some useful FMan Files */

#define FILE_SET_ENERGY       72    // Excitation energies CIS states
#define FILE_ENERGY           99    // 
#define FILE_NUCLEAR_GRADIENT 131   //
#define FILE_EFIELD           329   // 
#define FILE_DERCOUP          967   // Derrivative coupling 

/* And some offsets */
#define FILE_POS_CRNT_TOTAL_ENERGY  11
#define FILE_POS_BEGIN              0

#endif
