#ifndef GIFS_SH_QM_QCHEM_HPP
#define GIFS_SH_QM_QCHEM_HPP

#include "properties.hpp"
#include "qm_interface.hpp"
#include <vector>
#include <map>
#include <armadillo>

using REMKeys = std::map<std::string, std::string>;

class QM_QChem: public QMInterface{
public:
  QM_QChem(const std::vector<int> &qmid, int charge, int mult);
  void get_properties(PropMap &props);
  
private:  
  void get_nac_vector(arma::mat *nac, size_t A, size_t B);
  void get_wf_overlap(arma::mat *U);

  void get_gradient(arma::mat &g_qm, arma::uword surface);
  void get_gradient(arma::mat &g_qm);
  void get_gradient(arma::mat &g_qm, arma::mat &g_mm);
  void get_gradient(arma::mat &g_qm, arma::mat &g_mm, arma::uword surface);

  void get_ground_energy(arma::vec *e);
  void get_all_energies(arma::vec *e);

  std::ofstream get_input_handle(void);
  void write_molecule_section(std::ostream &ifile);
  void write_rem_section(std::ostream &os, const REMKeys &options);
  REMKeys excited_rem(void);
  void exec_qchem(void);

  void parse_qm_gradient(arma::mat &g_qm);
  void parse_mm_gradient(arma::mat &g_mm);
  void parse_energies(arma::vec &e);
  
  size_t readQFMan(int filenum, double * memptr, size_t N, size_t offset);
  
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

  int e_call_idx = -1;  // call index for ground and (e)xcited energy
  int ee_call_idx = -1;
};

/* Some useful FMan Files */

#define FILE_SET_ENERGY       72    // CIS state *excitation* energies 
#define FILE_ENERGY           99    // 
#define FILE_NUCLEAR_GRADIENT 131   //
#define FILE_EFIELD           329   // 
#define FILE_DERCOUP          967   // Derrivative coupling 
#define FILE_WF_OVERLAP       398   // wavefunction overlap

/* And some offsets */
#define FILE_POS_CRNT_TOTAL_ENERGY  11
#define FILE_POS_BEGIN              0

#endif
