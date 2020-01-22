#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>
#include <limits>
#include "qm_qchem.hpp"
#include "properties.hpp"

#define FMT "%12.8g"

QM_QChem::QM_QChem(const std::vector<int> &qmid, int charge, int mult):
  QMInterface(qmid, charge, mult),
  qc_scratch_directory(get_qcscratch()),
  qc_executable(get_qcprog()),
  exchange_method("HF"), // FIXME: Method/basis should be configurable;
  basis_set("6-31+G*"),  // by parsing an input file? Or higher up?
  excited_states(0)
{
  /*
    FIXME: Threading should be done in a configurable fashion
    set the number of threads, but don't overwrite if the flag is set elsewhere
  */
  setenv("QCTHREADS", "4", 0);

  /*
    FIXME: Setting $QCTHREADS seems sufficient on the subotnik
    cluster, but perhaps this is unique to our setup. Need to check
    with Evgeny to be sure.
  */
  //setenv("OMP_NUM_THREADS", "4", 0);
}


void QM_QChem::get_properties(PropMap &props){
  std::vector<double>& g_qm=props.get(QMProperty::qmgradient);
  std::vector<double>& g_mm=props.get(QMProperty::mmgradient);
  std::vector<double>& e=props.get(QMProperty::energies);

  (void) g_mm;
  (void) e;
  get_gradient_energies(g_qm, g_mm, e);
  //get_excited_gradient(g_qm, g_mm, e, 2);
  //get_nac_vector(g_qm, 1, 2);
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
  //FIXME: need to verify surface <= excited_states
  
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

  if (excited_states > 0){
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
    std::cout << e.size() << std::endl;
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
  size_t count = readQFMan(FILE_EFIELD, g_mm, 0, 3*NMM);

  /*
    FIXME: Why does the number of MM atoms fluctuate during a qmmm
    run?
  */
  if (count > 3 * NMM){
    g_mm.resize(3*NMM);
  }
  else if (count < 3* NMM){
    throw std::runtime_error("Unable to parse MM gradient!");
  }

  std::cout << std::endl;
  
  for (size_t i = 0; i < NMM; i++){
    const double q = chg_mm[i];
    g_mm[3*i + 0] *= q;
    g_mm[3*i + 1] *= q;
    g_mm[3*i + 2] *= q;
  }
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

// void QM_QChem::parse_mm_gradient(std::vector<double> &g_mm){
//   /*
//     efield.dat uses atomic units; its format is:

//     Ex(mm1) Ey(mm1) Ez(mm1)
//     ...
//     Ex(mmN) Ey(mmN) Ez(mmN)
//     Ex(qm1) Ey(qm1) Ez(qm1)
//     ...
//     Ex(qmN) Ey(qmN) Ez(qmN)

//     where Ea(u) is the component of the electric field in the a
//     direction at u and mmi and qmi are the coordinates of the ith mm
//     and qm atoms respectively.
//   */


//   /*
//     FIXME: Should be doing binary io with FMan for the efield.
//   */
  
//   /*
//     FIXME: This will return a *force* rather than a gradient unless it
//     is the negative of the field, which it may be. We need to return a
//     gradient.
//   */

//   std::string gfile = qc_scratch_directory + "/" + "efield.dat";
//   std::ifstream ifile(gfile);
//   std::string line;

//   size_t i = 0;
//   {
//     double x, y, z;
//     while (getline(ifile, line) && i < NMM){
//       std::stringstream s(line);
//       s >> x >> y >> z;
//       const double q = chg_mm[i];
//       g_mm[i*3 + 0] = x * q;
//       g_mm[i*3 + 1] = y * q;
//       g_mm[i*3 + 2] = z * q;
//       i++;
//     }
//   }

//   if (i != NMM){
//     throw std::runtime_error("Unable to parse MM gradient!");
//   }
// }



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

void QM_QChem::write_rem_section(std::ostream &os, std::map<std::string, std::string> options){
  // Default options 
  std::map<std::string, std::string> rem_keys
    {
     {"method",        exchange_method},
     {"basis",         basis_set},
     {"sym_ignore",    "true"},
     {"qm_mm",         "true"},
     {"qmmm_ext_gifs", "1"},
     {"input_bohr",    "true"}
    };

  // merge in options; usually will include jobtype 
  rem_keys.insert(options.begin(), options.end());
  
  if (! first_call ){
    rem_keys.emplace("scf_guess", "read");
  }  
  
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
