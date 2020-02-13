#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <limits>
#include <algorithm>
#include "qm_qchem.hpp"
#include "properties.hpp"
#include <armadillo>



QM_QChem::QM_QChem(const std::vector<int> &qmid, int charge, int mult):
  QMInterface(qmid, charge, mult),
  qc_scratch_directory(get_qcscratch()),
  qc_executable(get_qcprog()),
  exchange_method("HF"), // FIXME: Method/basis should be configurable;
  basis_set("6-31+G*"),  // by parsing an input file? Or higher up?
  excited_states(0)
{
  /*
    Set the number of threads, but don't overwrite if the flag is set
    elsewhere
  */
  // FIXME: Threading should be done in a configurable fashion
  setenv("QCTHREADS", "4", 0);
  /*
    Setting $QCTHREADS seems sufficient on the Subotnik cluster, but
    perhaps this is unique to our setup. Need to check with Evgeny to
    be sure. Changes to $OMP_NUM_THREADS seem to have no effect.
  */
}


void QM_QChem::get_properties(PropMap &props){
  std::cout << ":: [" << call_idx() << "]" << std::endl;
  // arma::cube U (excited_states, excited_states,1);
  // get_wf_overlap(&U);
  // U.slice(0).print();
  // std::cout << std::endl;

  // arma::cube D (3, (NMM + NQM),1);
  // get_nac_vector(&D, 1, 2);
  // D.slice(0).t().print();
  // std::cout << std::endl;

  for (QMProperty p: props.keys()){
    switch(p){

    case QMProperty::wfoverlap:
      get_wf_overlap(props.get(QMProperty::wfoverlap));
      break;

    case QMProperty::nacvector_imag:
      throw std::invalid_argument("Imaginary NAC not implemented!");
      break;

    case QMProperty::nacvector:{
      const arma::uvec &idx = *props.get_idx(QMProperty::nacvector);
      size_t A = idx[0]; size_t B = idx[1];
      get_nac_vector(props.get(QMProperty::nacvector), A, B);
      e_call_idx = call_idx(); ee_call_idx = call_idx();
      break;
    }

    case QMProperty::mmgradient:
      break; // if we need mm, then we will do qm, which catches both
    case QMProperty::qmgradient:{
      arma::cube *g_qm_c = props.get(QMProperty::qmgradient);
      arma::cube *g_mm_c = props.get(QMProperty::mmgradient);
      
      if (props.has_idx(QMProperty::qmgradient)){//excited states
        const arma::uvec &surf_idx = *props.get_idx(QMProperty::qmgradient);
        for (size_t i = 0; i < surf_idx.size(); i++){
	  arma::mat &g_qm = g_qm_c->slice(i);
	  if (props.has(QMProperty::mmgradient)){
	    arma::mat &g_mm = g_mm_c->slice(i);
	    get_excited_gradient(&g_qm, &g_mm, surf_idx[i]);
	  }
	  else{
	    get_excited_gradient(&g_qm, nullptr, surf_idx[i]);
	  }
        }
       	e_call_idx = call_idx(); ee_call_idx = call_idx();
      }
      else{ // ground states only
	arma::mat &g_qm = g_qm_c->slice(0);
	if (props.has(QMProperty::mmgradient)){
	  arma::mat &g_mm = g_mm_c->slice(0);
	  get_ground_gradient(&g_qm, &g_mm);
	}
	else{
	  get_ground_gradient(&g_qm, nullptr);
	}
	e_call_idx = call_idx();
      }
      break;
    }

    case QMProperty::energies:{
      // FIXME: clean-up the way we avoid an extra call for energy
      // (currently: sorting + e(e)_call_idx flags)
      if (props.has_idx(QMProperty::energies)){
      	get_all_energies(props.get(QMProperty::energies));
      }
      else{
      	get_ground_energy(props.get(QMProperty::energies)); // FIXME: combine ground/excited calls
      }
      break;
    }

    default: // Compile with the default case commented out and the compiler will detect neglected properties
      throw std::invalid_argument("Unknown QMProperty!");
      break;
    }
  }
}

void QM_QChem::get_wf_overlap(arma::cube *U){  
  static int last_overlap_call = -1;
  if (last_overlap_call != call_idx()){ //FIXME: make sure this is the right thing to do
    REMKeys k = excited_rem();
    k.insert({{"jobtype","sp"},
    	      {"namd_lowestsurface","1"},  //FIXME: parametrize lowest surface
    	      {"dump_wf_overlap", "1"}});

    std::ofstream input = get_input_handle();
    write_rem_section(input, k);
    write_molecule_section(input);
    input.close();
    exec_qchem();
  }
  else{
    //don't recompute 
  }
  
  last_overlap_call = call_idx();
  
  size_t count = readQFMan(FILE_WF_OVERLAP, U->memptr(), excited_states*excited_states, FILE_POS_BEGIN);
  if (count != excited_states * excited_states){
    throw std::runtime_error("Unable to parse wavefunction overlap");
  }  
}


void QM_QChem::get_ground_energy(arma::cube *e){
  // Build job if we need to
  if (e_call_idx != call_idx()){
    e_call_idx = call_idx();
    std::ofstream input = get_input_handle();
    write_rem_section(input, {{"jobtype","sp"}});
    write_molecule_section(input);
    input.close();
    exec_qchem();
  }

  parse_energies(e);
}

// Gets all energies: ground + excited states
void QM_QChem::get_all_energies(arma::cube *e){
  if (ee_call_idx != call_idx()){
    ee_call_idx = call_idx();
    REMKeys k = excited_rem();
    k.insert({{"jobtype","sp"}});
    
    std::ofstream input = get_input_handle();
    write_rem_section(input, k);
    write_molecule_section(input);
    input.close();
    exec_qchem();
  }

  parse_energies(e);
}


void QM_QChem::get_ground_gradient(arma::mat *g_qm, arma::mat *g_mm){
  // Build job
  std::ofstream input = get_input_handle();
  write_rem_section(input, {{"jobtype","force"}});
  write_molecule_section(input);
  input.close();
  
  exec_qchem();

  parse_qm_gradient(g_qm);
  if (g_mm && NMM > 0){
    parse_mm_gradient(g_mm);
  }
}

void QM_QChem::get_excited_gradient(arma::mat *g_qm, arma::mat *g_mm, arma::uword surface){
  if (surface > excited_states){
    throw std::invalid_argument("Requested surface not computed!");
  }

  REMKeys k = excited_rem();
  k.insert({{"jobtype","force"},
	    {"cis_state_deriv", std::to_string(surface)}});
  
  std::ofstream input = get_input_handle();
  write_rem_section(input, k);
  write_molecule_section(input);
  input.close();

  exec_qchem();

  (void) g_qm; (void) g_mm;
  parse_qm_gradient(g_qm);
  if (g_mm && NMM > 0){
    parse_mm_gradient(g_mm);
  }  
}


void QM_QChem::get_nac_vector(arma::cube *nac, size_t A, size_t B){

  REMKeys k = excited_rem();
  k.insert({{"jobtype","sp"},
	    {"calc_nac", "true"},
	    {"cis_der_numstate", "2"}});//Always, in our case, between 2 states
  
  std::ofstream input = get_input_handle();
  write_rem_section(input, k);

  /*
    Section format:
    $derivative_coupling
      comment line
      A B ...
    $end
    where A, B, ... are the states between which to compute the
    derivative coupling. The $rem section must include
    cis_der_numstates equal to the total number of states
    specified. The ground state is 0.

  */
  input << "$derivative_coupling" << std::endl;
  input << "comment" << std::endl;
  input << A << " " << B << std::endl;
  input << "$end" << std::endl;
  
  write_molecule_section(input);
  input.close();

  exec_qchem();

  size_t count = readQFMan(FILE_DERCOUP, nac->memptr(), 3 * (NMM + NQM), FILE_POS_BEGIN);
  if (count != 3 * (NMM + NQM)){
    throw std::runtime_error("Unable to parse NAC vector!");
  }
}



void QM_QChem::parse_energies(arma::cube *e_cube){
  //Bit of a nasty hack
  std::vector<double> e(e_cube->memptr(),e_cube->memptr() + e_cube->n_elem);
  
  readQFMan(FILE_ENERGY, e.data(), 1, FILE_POS_CRNT_TOTAL_ENERGY);

  if (excited_states > 0 && e.size() > 1){
    double e_ground = e[0];

    size_t count = readQFMan(FILE_SET_ENERGY, e.data(), excited_states, FILE_POS_BEGIN);
    if (count != excited_states){
      throw std::runtime_error("Unable to parse energies!");
    }
    
    /* Energies for higher states given in terms of exitations so we
       must add in the ground. We also need to shift the values so the
       ground state is at the begining. */

    // for (size_t i = excited_states - 1; i > 0; i--){
    //   e[i] = e[i-1] + e_ground;
    // }
    // e[0] = e_ground;

    e.pop_back();
    e.emplace(e.begin(), 0.0);
    for (auto & ei: e){
      ei += e_ground;
    }
  }
}


void QM_QChem::parse_qm_gradient(arma::mat *g_qm){
  size_t count = readQFMan(FILE_NUCLEAR_GRADIENT, g_qm->memptr(), 3*NQM, FILE_POS_BEGIN);
  if (count != 3 * NQM){
    throw std::runtime_error("Unable to parse QM gradient!");
  }
}


/*
  FIXME: This will return a *force* rather than a gradient unless
  efield is the negative of the field, which it may be. We need to
  return a gradient.
*/
void QM_QChem::parse_mm_gradient(arma::mat *g_mm){
  /* Set a ceiling on the number of doubles we read-in because the
     number of MM atoms can fluctuate during a simulation. */
  size_t count = readQFMan(FILE_EFIELD, g_mm->memptr(), 3*NMM, FILE_POS_BEGIN);

  if (count != 3 * NMM){
    throw std::runtime_error("Unable to parse MM gradient!");
  }

  for (size_t i = 0; i < NMM; i++){
    const double q = chg_mm[i];
    arma::vec f = g_mm->unsafe_col(i);
    f *= q;
  }
  
  // for (size_t i = 0; i < NMM; i++){
  //   const double q = chg_mm[i];
  //   g_mm[3*i + 0] *= q;
  //   g_mm[3*i + 1] *= q;
  //   g_mm[3*i + 2] *= q;
  // }
  /*
    This method replaced parsing efield.dat, which has the following
    format:

    Ex(mm1) Ey(mm1) Ez(mm1)
    ...
    Ex(mmN) Ey(mmN) Ez(mmN)
    Ex(qm1) Ey(qm1) Ez(qm1)
    ...
    Ex(qmN) Ey(qmN) Ez(qmN)

    where Ea(u) is the component of the electric field in the a
    direction at u and mmi and qmi are the coordinates of the ith mm and
    qm atoms respectively and there is a perhaps-present factor of (-1)
    s.t. the product of the value and the charge is a gradient rather
    than a force.
  */
}


const std::string QM_QChem::get_qcprog(void){
  char * qc_str = std::getenv("QC");
  std::string default_path = "qcprog.exe";
  if (nullptr == qc_str){
    std::cerr << "Warning, $QC not set; using " << default_path << std::endl;
    return default_path;
  }
  else{
    return std::string(qc_str) + "/exe/qcprog.exe";
  }
}


const std::string QM_QChem::get_qcscratch(void){
  char * pwd = std::getenv("PWD");
  std::string scratch_path = std::string(pwd) + "/GQSH.sav";
  
  char * qc_str = std::getenv("QCSCRATCH");
  if (nullptr == qc_str){
    std::cerr << "Warning, $QCSCRATCH not set; using " << scratch_path << std::endl;
  }
  else if (qc_str[0] != '/'){
    /*
      FIXME: *Q-Chem* will crash if QCSCRATCH=""
    */
    std::cerr << "Warning, $QCSCRATCH is not an absolute path; instead using " << scratch_path << std::endl;
  }
  else{
    scratch_path = std::string(qc_str);
  }

  /* Make sure the directory exists; create otherwise */
  struct stat st = {};
  if (-1 == stat(scratch_path.c_str(), &st)){
    mkdir(scratch_path.c_str(), 0700);
  }

  return scratch_path;
}


void QM_QChem::exec_qchem(void){
  std::string cmd = "cd " + qc_scratch_directory + "; " + // WD->$QCSCRATCH
    qc_executable + " " + qc_input_file + " " + qc_scratch_directory + " >" + qc_log_file;
  int status = std::system(cmd.c_str());
  if (status){
    throw std::runtime_error("Q-Chem could not be called or exited abnormally; see " + qc_log_file);
  }
  first_call = false;
}


/* Consistent formatting and file name */
std::ofstream QM_QChem::get_input_handle(){
  std::ofstream os(qc_scratch_directory + "/" + qc_input_file);
  os.setf(std::ios_base::fixed, std::ios_base::floatfield);
  os.precision(std::numeric_limits<double>::digits10);

  return os;
}

REMKeys QM_QChem::excited_rem(void){
  if (! (excited_states > 0)){
    throw std::logic_error("Cannot request multi-state property with out excited states.");
  }

  REMKeys excited
    {
     {"cis_n_roots", std::to_string(excited_states)},
     {"cis_singlets", "true"},  //Fixme: singlet/triplet selection should
     {"cis_triplets", "false"}  //be configurable
    };

  return excited;
}

void QM_QChem::write_rem_section(std::ostream &os, const REMKeys &options){
  // Default options 
  REMKeys rem_keys
    {
     {"method",        exchange_method},
     {"basis",         basis_set},
     {"sym_ignore",    "true"},
     {"qm_mm",         "true"},
     {"input_bohr",    "true"}
    };
  
  if (first_call ){
    rem_keys.emplace("qmmm_ext_gifs", "1");
  }
  else{
    rem_keys.emplace("qmmm_ext_gifs", "2");
    rem_keys.emplace("scf_guess",     "read");
  }  

  /*
    FIXME: make sure passed options take precedence?  Could do this by
    inserting the defaults into the passed options rather than the
    other way around.
  */
  
  // merge in options; usually will include jobtype 
  rem_keys.insert(options.begin(), options.end());
  
  os << "$rem" << std::endl;
  for (const auto& e: rem_keys){
    os << e.first << " " << e.second << std::endl;
  }
  os << "$end" << std::endl;  
}


void QM_QChem::write_molecule_section(std::ostream &os){
  /*
    format of $molecule & $externall_charges sections:
    $molecule
    [charge] [multiplicity]
    [atomic-number] [x-coord] [y-coord] [x-coord]
    ...
    $end

    $external_charges
    [x-coord] [y-coord] [x-coord] [charge]
    ...
    $end
  */  
  os << "$molecule" << std::endl;
  os << qm_charge << " " << qm_multiplicity << std::endl;

  for (size_t i = 0; i < NQM; i++){
    os << atomids[i]      << " ";        // id
    os << crd_qm[i*3 + 0] << " ";        // x
    os << crd_qm[i*3 + 1] << " ";        // y
    os << crd_qm[i*3 + 2] << std::endl;  // x
  }
  os <<  "$end" << std::endl;

  if (NMM > 0){
    os << std::endl << "$external_charges" << std::endl;    
    for (size_t i = 0; i < NMM; i++){
      os << crd_mm[i*3 + 0] << " ";        // x
      os << crd_mm[i*3 + 1] << " ";        // y
      os << crd_mm[i*3 + 2] << " ";        // z
      os << chg_mm[i]       << std::endl;  // charge
    }
    os << "$end" << std::endl;
  }
}


/*
  Given a q-qchem file number (see qm_qchem.hpp for examples), read N
  elements from the file (in the scratch directory) starting from the
  offset into memptr, which must point to a block of memory with at
  least N elements.
*/
size_t QM_QChem::readQFMan(int filenum, double * memptr, size_t N, size_t offset){  
  std::string path = qc_scratch_directory + "/" + std::to_string(filenum) + ".0";
  std::ifstream ifile;
  ifile.open(path, std::ios::in | std::ios::binary);

  if (!ifile.is_open()){
    throw std::ios_base::failure("Error: cannot read from " + path);
  }

  char buffer[sizeof(double)];
  size_t i = 0;
  size_t count = 0;
  while(ifile.read(buffer, sizeof(double))){
    if ((i >= offset) && (i < offset + N)){
      double * d = (double *) buffer;
      *(memptr + count++) = *d;
    }
    i++;
  }
  ifile.close();

  return count;
}
