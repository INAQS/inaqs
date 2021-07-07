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
           const arma::uvec& in_qmids, 
	   arma::mat& in_qm_crd, 
	   arma::mat& in_mm_crd, 
	   arma::vec& in_mm_chg, 
	   const int charge, 
	   const int mult,
	   const int excited_states,
           const int min_state);
  
  void get_properties(PropMap &props) override;

private:
  void state_tracker(PropMap &props);

  void get_nac_vector(arma::mat *nac, size_t A, size_t B);
  void get_wf_overlap(arma::mat *U);
  void get_diabatic_rot_mat(arma::mat *U);

  void get_gradient(arma::mat &g_qm, arma::uword surface);
  void get_gradient(arma::mat &g_qm, arma::mat &g_mm, arma::uword surface);

  void get_ground_energy(arma::vec *e);
  void get_all_energies(arma::vec *e);

  void do_state_analysis(void);
  void do_record_spectrum(void);
  void do_boys_diabatization(void);
  
  std::ofstream get_input_handle(void);
  void write_molecule_section(std::ostream &ifile);
  void write_rem_section(std::ostream &os, const REMKeys &options);
  REMKeys excited_rem(void);
  void exec_qchem(void);

  void parse_qm_gradient(arma::mat &g_qm);
  void parse_mm_gradient(arma::mat &g_mm);
  void parse_energies(arma::vec &e);

  
  
  size_t readQFMan(int filenum, double * memptr, size_t N, size_t offset);
  template<typename T>
  void readQFMan(int filenum, T & a, size_t offset=0){
    if (a.n_elem != readQFMan(filenum, a.memptr(), a.n_elem, offset)){
      throw std::logic_error("Failure to parse file " + std::to_string(filenum));
    }
  }

  ConfigBlockReader qchem_reader();
  
  const std::string get_qcprog(void);
  const std::string get_qcscratch(std::string conf_dir);
  const std::string get_qcwdir(void);


  // FIXME: should refactor these into an options class
  std::string qc_scratch_directory;
  std::string qc_executable;
  std::string qc_input_file = "GQSH.in";
  std::string qc_log_file = "QCHEM.out";
  std::string exchange_method;
  std::string scf_algorithm;
  std::string basis_set;
 
  bool first_call = true;

  bool singlets = true;  // Defaults for CIS calculation
  bool triplets = false;
  bool spin_flip = false;
  bool track_states = false;

  size_t shstates = 0;
  
  bool state_analysis = false;
  bool save_nacvector = false;
  bool record_spectrum = false;
  std::vector<int> boys_states {};
  bool boys_diabatization = false;
  
  
  enum class S{
    energy,
    ex_energy,
    ex_grad,
    wfoverlap,
    once,
  };

  // As long as QM_Interface.update() is called before
  // get_properites(), call_idx() will always return a value > 0. This
  // is required behavior.
  bool called(S s);
};

/* Some useful FMan Files */

#define FILE_SET_ENERGY        72    // CIS state *excitation* energies
#define FILE_ENERGY            99    //
#define FILE_NUCLEAR_GRADIENT 131    //
#define FILE_EFIELD           329    //
#define FILE_DERCOUP          967    // Derrivative coupling
#define FILE_WF_OVERLAP       398    // wavefunction overlap
#define FILE_DIAB_ROT_MAT     941    // for diabatization rotations
#define FILE_TRANS_DIP_MOM    942    // Transition dipole moments: states along cols, rows: strength, x, y, z
#define FILE_CIS_S2          1200

/* And some offsets */
#define FILE_POS_CRNT_TOTAL_ENERGY  11
#define FILE_POS_BEGIN              0

#endif
