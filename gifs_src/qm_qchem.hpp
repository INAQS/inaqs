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
	   const int excited_states_in,
           const int min_state_in);
  
  void get_properties(PropMap &props) override;

private:
  enum class QCFILE{
    FILE_SET_ENERGY       =   72,    // CIS state *excitation* energies
    FILE_ENERGY           =   99,    //
    FILE_NUCLEAR_GRADIENT =  131,    //
    FILE_EFIELD           =  329,    //
    FILE_DERCOUP          =  967,    // Derrivative coupling
    FILE_WF_OVERLAP       =  398,    // wavefunction overlap
    FILE_DIAB_ROT_MAT     =  941,    // for diabatization rotations
    FILE_TRANS_DIP_MOM    =  942,    // Transition dipole moments: states along cols, rows: strength, x, y, z
    FILE_CIS_S2           = 1200,    // S^2 matrix elements for CIS states
  };
  std::string to_string(const QCFILE f);

  void state_tracker(PropMap &props);

  void get_nac_vector(arma::mat &nac, size_t A, size_t B);
  void get_wf_overlap(arma::mat &U);
  void get_diabatic_rot_mat(arma::mat &U);

  void get_gradient(arma::mat &g_qm, arma::uword surface);
  void get_gradient(arma::mat &g_qm, arma::mat &g_mm, arma::uword surface);

  void get_ground_energy(arma::vec &e);
  void get_all_energies(arma::vec &e);

  void do_state_analysis(void);
  void do_record_spectrum(void);
  void do_boys_diabatization(void);
  
  std::ofstream get_input_handle(void);
  void write_molecule_section(std::ostream &ifile);
  void write_rem_section(std::ostream &os, const REMKeys &options);
  REMKeys excited_rem(void);
  void exec_qchem(void);
  std::string qc_log_file_idx(void);

  void parse_mm_gradient(arma::mat &g_mm);
  void parse_energies(arma::vec &e);
  
  size_t readQFMan(QCFILE filenum, double * memptr, size_t N, size_t offset);
  template<typename T> // template for arma tensors
  void readQFMan(QCFILE filenum, T & a, size_t offset=0){
    if (a.n_elem != readQFMan(filenum, a.memptr(), a.n_elem, offset)){
      throw std::logic_error("Failure to parse file: " + to_string(filenum));
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
  bool dump_qc_output = false;
  bool dump_qc_input = false;

  double sing_thresh = 0;
  
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

  // Some offsets for files
  const size_t POS_BEGIN              = 0;
  const size_t POS_CRNT_TOTAL_ENERGY  = 11;
};

#endif
