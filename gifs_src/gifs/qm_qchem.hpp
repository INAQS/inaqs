#ifndef GIFS_SH_QM_QCHEM_HPP
#define GIFS_SH_QM_QCHEM_HPP

#include "properties.hpp"
#include "qm_interface.hpp"
#include "configreader.hpp"
#include <vector>
#include <map>
#include <armadillo>

using REMKeys = std::map<std::string, std::string>;

class QM_QChem: public QMInterface{
public:
  QM_QChem(FileHandle& fh, 
       arma::uvec& in_qmids, 
	   arma::mat& in_qm_crd, 
	   arma::mat& in_mm_crd, 
	   arma::vec& in_mm_chg, 
	   int charge, 
	   int mult,
	   size_t excited_states);
  void get_properties(PropMap &props);

private:  
  void get_nac_vector(arma::mat *nac, size_t A, size_t B);
  void get_wf_overlap(arma::mat *U);

  void get_gradient(arma::mat &g_qm, arma::uword surface);
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

  ConfigBlockReader qchem_reader();
  
  const std::string get_qcprog(void);
  const std::string get_qcscratch(std::string conf_dir);

  std::string qc_scratch_directory;
  std::string qc_executable;
  std::string qc_input_file = "GQSH.in";
  std::string qc_log_file = "GQSH.out";
  std::string exchange_method;
  std::string basis_set;
  bool first_call = true;

  enum class S{
    energy,
    ex_energy,
    ex_grad,
    wfoverlap,
  };

  //FIXME: when switching to the new interface, verify that call_idx()
  //will always return a value > 0
  bool called(S s);
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
