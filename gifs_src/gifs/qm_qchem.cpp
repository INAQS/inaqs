#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <limits>
#include <algorithm>
#include "qm_qchem.hpp"
#include "properties.hpp"



QM_QChem::QM_QChem(const std::vector<int> &qmid, int charge, int mult):
  QMInterface(qmid, charge, mult),
  qc_scratch_directory(get_qcscratch()),
  qc_executable(get_qcprog()),
  exchange_method("HF"), // FIXME: Method/basis should be configurable;
  basis_set("6-31+G*"),  // by parsing an input file? Or higher up?
  excited_states(0)
{
  /*
    FIXME: Threading should be done in a configurable fashion set the
    number of threads, but don't overwrite if the flag is set
    elsewhere
  */
  setenv("QCTHREADS", "4", 0);

  /*
    Setting $QCTHREADS seems sufficient on the Subotnik cluster, but
    perhaps this is unique to our setup. Need to check with Evgeny to
    be sure. Changes to $OMP_NUM_THREADS seem to have no effect.
  */
}

void QM_QChem::get_properties(PropMap &props){
  const std::vector<int> *idx;
  
  for (QMProperty p: props.keys()){
    switch(p){
      
    case QMProperty::wfoverlap:
      get_wf_overlap(*props.get(QMProperty::wfoverlap));
      break;
      
    case QMProperty::nacvector_imag:
      throw std::invalid_argument("Imaginary NAC not implemented!");
      break;
      
    case QMProperty::nacvector:{
      idx = props.get_idx(QMProperty::nacvector);
      size_t A = (*idx)[0]; size_t B = (*idx)[1];
      get_nac_vector(*props.get(QMProperty::nacvector), A, B);
      e_call_idx = call_idx(); ee_call_idx = call_idx();
      break;
    }
      
    case QMProperty::mmgradient:
      break; // if we need mm, then we will do qm, which catches both
    case QMProperty::qmgradient:
      if (props.has_idx(QMProperty::qmgradient)){
	const std::vector<int> &surf_idx = *props.get_idx(QMProperty::qmgradient);
	std::vector<double> &qm_gradients = *props.get(QMProperty::qmgradient);
	std::vector<double> &mm_gradients = *props.get(QMProperty::mmgradient);
	for (size_t i = 0; i < surf_idx.size(); i++){ // FIXME: would like to do this without copying; do vector "views" exist?
	  std::vector<double> qmg(qm_gradients.begin() + i*3*NQM, qm_gradients.begin() + (i+1)*3*NQM);
	  std::vector<double> mmg(mm_gradients.begin() + i*3*NMM, mm_gradients.begin() + (i+1)*3*NMM);
	  
	  get_excited_gradient(&qmg, &mmg, (size_t) surf_idx[i]);
	  
	  std::copy(qmg.begin(), qmg.end(), qm_gradients.begin() + i*3*NQM);
	  std::copy(mmg.begin(), mmg.end(), mm_gradients.begin() + i*3*NMM);
	}
	e_call_idx = call_idx(); ee_call_idx = call_idx();
      }
      else{
	get_ground_gradient(props.get(QMProperty::qmgradient), // FIXME: combine ground/excited calls
			    props.get(QMProperty::mmgradient));
	e_call_idx = call_idx();
      }
      break;

    case QMProperty::energies:
      // FIXME: clean-up the way we avoid an extra call for energy
      // (currently: sorting + e(e)_call_idx flags)
      if (props.has_idx(QMProperty::energies)){
	get_all_energies(props.get(QMProperty::energies));
      }
      else{
	get_ground_energy(props.get(QMProperty::energies)); // FIXME: combine ground/excited calls
      }
      break;

    default: // Compile with the default case commented out and the compiler will detect neglected properties
      throw std::invalid_argument("Unknown QMProperty!");
      break;
    }
  }
}

//FIXME: make sure this is only computed on demand => modify q-chem
void QM_QChem::get_wf_overlap(std::vector<double> &U){ //FIXME: use an armadillo matrix
  static int overlap_call = 0;
  if (overlap_call != call_idx()){ //FIXME: make sure this is the right thing to do
    overlap_call = call_idx();
    std::ofstream input = get_input_handle();
    //FIXME: should compose CIS functionality rather than duplicating
    write_rem_section(input,
		      {{"jobtype","sp"},
		       {"namd_lowestsurface","1"},  //FIXME: parametrize lowest surface
		       {"dump_wf_overlap", "1"},
		       {"cis_n_roots", std::to_string(excited_states)},
		       {"cis_singlets", "true"},    //Fixme: singlet/triplet selection should
		       {"cis_triplets", "false"}}); //be configurable
    write_molecule_section(input);
    input.close();
    exec_qchem();
  }
  else{
    overlap_call = call_idx();
    //don't recompute 
  }
  
  size_t count = readQFMan(FILE_WF_OVERLAP, U);
  if (count != excited_states * excited_states){
    throw std::runtime_error("Unable to parse wavefunction overlap");
  }
  
}


void QM_QChem::get_ground_energy(std::vector<double> *e){
  // Build job if we need to
  if (e_call_idx != call_idx()){
    std::ofstream input = get_input_handle();
    write_rem_section(input, {{"jobtype","sp"}});
    write_molecule_section(input);
    input.close();
    exec_qchem();
  }
  
  parse_energies(*e);
}

// Gets all energies: ground + excited states
void QM_QChem::get_all_energies(std::vector<double> *e){
  if (ee_call_idx != call_idx()){
    std::ofstream input = get_input_handle();
    write_rem_section(input,
		      {{"jobtype","sp"},
		       {"cis_n_roots", std::to_string(excited_states)},
		       {"cis_singlets", "true"},  //Fixme: singlet/triplet selection should
		       {"cis_triplets", "false"}  //be configurable
		      });
    write_molecule_section(input);
    input.close();
    exec_qchem();
  }

  parse_energies(*e);
}


void QM_QChem::get_ground_gradient(std::vector<double> *g_qm,
				   std::vector<double> *g_mm){
  // Build job
  std::ofstream input = get_input_handle();
  write_rem_section(input, {{"jobtype","force"}});
  write_molecule_section(input);
  input.close();
  
  exec_qchem();

  parse_qm_gradient(*g_qm);
  if (g_mm && NMM > 0){
    parse_mm_gradient(*g_mm);
  }
}

void QM_QChem::get_excited_gradient(std::vector<double> *g_qm,
				    std::vector<double> *g_mm,
				    size_t surface){
  if (surface > excited_states){
    throw std::invalid_argument("Requested surface not computed!");
  }
  
  std::ofstream input = get_input_handle();
  write_rem_section(input,
		    {{"jobtype","force"},
		     {"cis_n_roots", std::to_string(excited_states)},
		     {"cis_state_deriv", std::to_string(surface)},
		     {"cis_singlets", "true"},  //Fixme: singlet/triplet selection should
		     {"cis_triplets", "false"}  //be configurable
		    });
  write_molecule_section(input);
  input.close();

  exec_qchem();

  parse_qm_gradient(*g_qm);
  if (g_mm && NMM > 0){
    parse_mm_gradient(*g_mm);
  }  
}


void QM_QChem::get_gradient_energies(std::vector<double> &g_qm,
				     std::vector<double> &g_mm,
				     std::vector<double> &e){

  // Build job
  std::ofstream input = get_input_handle();
  write_rem_section(input, {{"jobtype","force"}});
  write_molecule_section(input);
  input.close();
  
  exec_qchem();

  parse_energies(e);
  parse_qm_gradient(g_qm);
  if (NMM > 0){
    parse_mm_gradient(g_mm);
  }
}


void QM_QChem::get_excited_gradient(std::vector<double> &g_qm,
				    std::vector<double> &g_mm,
				    std::vector<double> &e, size_t surface){
  if (surface > excited_states){
    throw std::invalid_argument("Requested surface is greater than the number of excited states!");
  }
  
  std::ofstream input = get_input_handle();
  write_rem_section(input,
		    {{"jobtype","force"},
		     {"cis_n_roots", std::to_string(excited_states)},
		     {"cis_state_deriv", std::to_string(surface)},
		     {"cis_singlets", "true"},  //Fixme: singlet/triplet selection should
		     {"cis_triplets", "false"}  //be configurable
		    });
  write_molecule_section(input);
  input.close();

  exec_qchem();

  parse_energies(e);
  parse_qm_gradient(g_qm);
  if (NMM > 0){
    parse_mm_gradient(g_mm);
  }  
}


void QM_QChem::get_nac_vector(std::vector<double> &nac, size_t A, size_t B){
  std::ofstream input = get_input_handle();
  write_rem_section(input,
		    {{"jobtype","sp"},
		     {"calc_nac", "true"},
		     {"cis_n_roots", std::to_string(excited_states)},
		     {"cis_der_numstate", "2"}, //Always, in our case, between 2 states
		     {"cis_singlets", "true"},  //Fixme: singlet/triplet selection should
		     {"cis_triplets", "false"}  //be configurable
		    });

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

  parse_nac_vector(nac);  
}


void QM_QChem::parse_nac_vector(std::vector<double> &nac){
  size_t count = readQFMan(FILE_DERCOUP, nac);
  if (count != 3 * (NMM + NQM)){
    throw std::runtime_error("Unable to parse NAC vector!");
  }
}


void QM_QChem::parse_energies(std::vector<double> &e){
  readQFMan(FILE_ENERGY, e, 1, FILE_POS_CRNT_TOTAL_ENERGY);

  if (excited_states > 0 && e.size() > 1){
    double e_ground = e[0];

    readQFMan(FILE_SET_ENERGY, e);
    /* Energies for higher states given in terms of exitations so we
       must add in the ground */
    e.emplace(e.begin(), 0.0);
    for (auto & ei: e){
      ei += e_ground;
    }
  }

  if (e.size() != excited_states + 1){
    throw std::runtime_error("Unable to parse energies!");
  }
}


void QM_QChem::parse_qm_gradient(std::vector<double> &g_qm){
  size_t count = readQFMan(FILE_NUCLEAR_GRADIENT, g_qm);
  if (count != 3 * NQM){
    throw std::runtime_error("Unable to parse QM gradient!");
  }
}


/*
  FIXME: This will return a *force* rather than a gradient unless
  efield is the negative of the field, which it may be. We need to
  return a gradient.
*/
void QM_QChem::parse_mm_gradient(std::vector<double> &g_mm){
  /* Set a ceiling on the number of doubles we read-in because the
     number of MM atoms can fluctuate during a simulation. */
  readQFMan(FILE_EFIELD, g_mm, 3*NMM, FILE_POS_BEGIN);

  if (g_mm.size() != 3 * NMM){
    throw std::runtime_error("Unable to parse MM gradient!");
  }

  for (size_t i = 0; i < NMM; i++){
    const double q = chg_mm[i];
    g_mm[3*i + 0] *= q;
    g_mm[3*i + 1] *= q;
    g_mm[3*i + 2] *= q;
  }
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
  std::string cmd = "cd " + qc_scratch_directory + "; " +
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


void QM_QChem::write_rem_section(std::ostream &os, const std::map<std::string, std::string> &options){
  // Default options 
  std::map<std::string, std::string> rem_keys
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
  Given a q-qchem file number (see qm_qchem.hpp for examples), reads
  the file from the scratch directory into the vector v, the contents
  of which are destroyed.
*/
size_t QM_QChem::readQFMan(int filenum, std::vector<double> &v){
  return readQFMan(filenum, v, v.max_size(), 0);
}

//FIXME: readQFMan should also be able to take an armadillo matrix/whatever/...

/*
  Same as above except we read N elements starting from offset.
*/
size_t QM_QChem::readQFMan(int filenum, std::vector<double> &v, size_t N, size_t offset){
  std::string path = qc_scratch_directory + "/" + std::to_string(filenum) + ".0";
  std::ifstream ifile;
  ifile.open(path, std::ios::in | std::ios::binary);

  if (!ifile.is_open()){
    throw std::ios_base::failure("Error: cannot read from " + path);
  }

  /*
    clear vector before reading so we aren't polluted by previous
    elements; probably a small performance hit here, but 1) it doesn't
    matter (this isn't our bottle-neck) and 2) it makes the code a lot
    clearer.
  */
  v.clear();
  char buffer[sizeof(double)];
  size_t i = 0;
  while(ifile.read(buffer, sizeof(double))){
    double * d = (double *) buffer;
    if ((i >= offset) && (i < offset+N)){
      v.push_back(*d);
    }
    i++;
  }
  ifile.close();

  return i;
}
