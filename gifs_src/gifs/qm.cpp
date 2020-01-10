#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
#include "qm.hpp"
#include "properties.hpp"

#define FMT "%12.8g"

QMInterface::QMInterface(size_t nqm, std::vector<int> &qmid):
  qc_scratch_directory(get_qcscratch()),
  qc_executable(get_qcprog()){
  NQM = nqm;
  NMM = 0;
  atomids = qmid;
  crd_qm.reserve(3 * nqm);

  qm_charge = 0;
  qm_multiplicity = 1;

  /*
    FIXME: This should be done in a configurable fashion
    set the number of threads, but don't overwrite if the flag is set elsewhere
  */
  setenv("QCTHREADS", "4", 0);
}

void QMInterface::update(std::vector<double> &crdqm,
			 std::vector<double> &crdmm,
			 std::vector<double> &chgmm){
  /*
    FIXME: the assignment operator '=' uses copying; when we implment
    our whole class hierarcy, memory will be owned one level above the
    BOMD class.
  */
  crd_qm = crdqm;
  NMM = chgmm.size();
  crd_mm = crdmm;
  chg_mm = chgmm;
}

void QMInterface::update(const float* crdqm, size_t nmm, const float* crdmm, const float* chgmm)
{
  NMM = nmm;
  if (chg_mm.capacity() < nmm){
    chg_mm.resize(nmm);
    crd_mm.resize(3 * nmm);
  }
  // QM Crd
  for (size_t i=0; i<NQM; ++i) {
      crd_qm[i*3    ] = crdqm[i*3    ] * 10;
      crd_qm[i*3 + 1] = crdqm[i*3 + 1] * 10;
      crd_qm[i*3 + 2] = crdqm[i*3 + 2] * 10;
  }
  // MM Crd
  for (size_t i=0; i<NMM; ++i) {
      crd_mm[i*3    ] = crdmm[i*3    ] * 10;
      crd_mm[i*3 + 1] = crdmm[i*3 + 1] * 10;
      crd_mm[i*3 + 2] = crdmm[i*3 + 2] * 10;
  }
  // Charges
  for (size_t i=0; i<NMM; ++i) {
      chg_mm[i] = chgmm[i];
  }
}

void QMInterface::get_properties(PropMap &props){
  std::vector<double>& g_qm=props.get(QMProperty::qmgradient);
  std::vector<double>& g_mm=props.get(QMProperty::mmgradient);
  std::vector<double>& e=props.get(QMProperty::energies);

  //auto& g_qm = *props[QMProperty::qmgradient];
  //auto& g_mm = *props[QMProperty::mmgradient];
  
  get_gradient_energies(g_qm, g_mm, e);
}

void QMInterface::get_gradient_energies(std::vector<double> &g_qm,
					std::vector<double> &g_mm,
					std::vector<double> &e){
  std::string ifname = "GQSH.in";
  std::string qcprog = get_qcprog();
  std::string savdir = "./GQSH.sav";

  write_gradient_job();
  exec_qchem();
  parse_qm_gradient(g_qm, e);
  if (NMM > 0){
    parse_mm_gradient(g_mm);
  }
}

std::string QMInterface::get_qcprog(void){
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

/*
  FIXME: need to ensure $QCSCRATCH is an absolute path
*/
std::string QMInterface::get_qcscratch(void){
  char * pwd = std::getenv("PWD");
  std::string default_path = std::string(pwd) + "/GQSH.sav";

  return default_path;
  
  char * qc_str = std::getenv("QCSCRATCH");
  if (nullptr == qc_str){
    std::cerr << "Warning, $QCSCRATCH not set; using " << default_path << std::endl;
    return default_path;
  }
  else{
    return std::string(qc_str);
  }
}

void QMInterface::parse_mm_gradient(std::vector<double> &g_mm){
  /*
    FIXME: need to verify the units of efield.dat. Here we assume they
    are: Hartrees/Angstrom/e-
  */

  std::string gfile = qc_scratch_directory + "/" + "efield.dat";
  std::ifstream ifile(gfile);
  std::string line;

  size_t i = 0;
  {
    double x, y, z;
    while (getline(ifile, line) && i < NMM){
      std::stringstream s(line);
      s >> x >> y >> z;
      const double q = chg_mm[i];
      g_mm[i*3 + 0] = x * q;
      g_mm[i*3 + 1] = y * q;
      g_mm[i*3 + 2] = z * q;
      i++;
    }
  }

  if (i != NMM){
    throw std::runtime_error("Unable to parse MM gradient!");
  }

  ang2bohr(g_mm);
}

void QMInterface::parse_qm_gradient(std::vector<double> &g_qm,
				    std::vector<double> &e){
  std::string gfile = qc_scratch_directory + "/" + "GRAD";
  std::ifstream ifile(gfile);
  std::string line;

  bool gradient = false;
  bool energy = false;
  size_t i = 0;
  
  while(getline(ifile, line) && i < NQM){
    std::stringstream s(line);
    if (energy){
      s>>e[0];
      energy = false;
    }

    if (gradient){
      s>>g_qm[i*3 + 0]
       >>g_qm[i*3 + 1]
       >>g_qm[i*3 + 2];
      i++;
    }
    
    if (line == "$energy"){
      energy = true;
    }
    if (line == "$gradient"){
      energy = false;
      gradient = true;
    }
  }

  if (i != NQM){
    throw std::runtime_error("Unable to parse QM gradient!");
  }
  
  ang2bohr(g_qm);
}


/*
  Unit conversion: Q-Chem outputs (??? need to confirm) in
  Hartree/Angstrom; our interface expects Hartree/Bohr.
*/
inline void QMInterface::ang2bohr(std::vector<double> &v){
  const double a2b = 0.529177249;
  for (auto& e : v){
    e *= a2b;
  }
}

void QMInterface::exec_qchem(void){
  std::string cmd = "cd " + qc_scratch_directory + "; " + qc_executable + " " + qc_input_file + " " + qc_scratch_directory;
  std::cout << cmd << std::endl;
  int status = std::system(cmd.c_str());
  if (status){
    throw std::runtime_error("Q-Chem could not be called or exited abnormally.");
  }
  first_call = false;
}

void QMInterface::write_gradient_job(){
  std::ofstream ifile(qc_scratch_directory + "/" + qc_input_file);
  ifile.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
  ifile << R"(
$rem
  jobtype          force
  method           hf
  basis	           6-31+G*
  sym_ignore       true
  qm_mm            true     # external charges in NAC; generate efield.dat
  qmmm_ext_gifs    1
)";

  /*
    NOTE:
    can add
      input_bohr = true
    to the q-chem input if our input is in atomic units.

    The default in angstroms and it's not clear how all of the
    gradient methods interact with a input unit change...
   */
  
  if (! first_call ){
    ifile << "  scf_guess        read" << std::endl;
  }
  
  ifile << "$end" << std::endl;

  /*
    format of $molecule section:
    $molecule
    [charge] [multiplicity]
    [atomic-number] [x-coord] [y-coord] [x-coord]
    ...
    $end
  */
  ifile << "$molecule" << std::endl;
  ifile << qm_charge << " " << qm_multiplicity << std::endl;
  /*
    FIXME: need to increase the precision of these outputs
  */
  for (size_t i = 0; i < NQM; i++){
    ifile << atomids[i]      << " ";        // id
    ifile << crd_qm[i*3 + 0] << " ";        // x
    ifile << crd_qm[i*3 + 1] << " ";        // y
    ifile << crd_qm[i*3 + 2] << std::endl;  // x
  }
  ifile <<  "$end" << std::endl;

  if (NMM > 0){
    ifile << std::endl << "$external_charges" << std::endl;
    /*
      format of $external_charges section:
      [x-coord] [y-coord] [x-coord] [charge]
    */
    
    for (size_t i = 0; i < NMM; i++){
      ifile << crd_mm[i*3 + 0] << " ";        // x
      ifile << crd_mm[i*3 + 1] << " ";        // y
      ifile << crd_mm[i*3 + 2] << " ";        // z
      ifile << chg_mm[i]       << std::endl;  // charge
    }
    ifile << "$end" << std::endl;
  }

  ifile.close();
}


int QMInterface::readQFMan(int filenum, std::vector<double> &v){
  std::string path = qc_scratch_directory + "/" + std::to_string(filenum) + ".0";
  std::ifstream ifile;
  ifile.open(path, std::ios::in | std::ios::binary);

  if (!ifile.is_open()){
    throw std::ios_base::failure("Error: cannot read from " + path);
  }
  
  char buffer[sizeof(double)];
  size_t i = 0;
  while(ifile.read(buffer, sizeof(double))){
    double * d = (double *) buffer;
    v.push_back(*d);
    i++;
  }
  ifile.close();

  return i;
}
