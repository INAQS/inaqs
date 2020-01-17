#include <vector>
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
  qc_executable(get_qcprog()){

  /*
    FIXME: This should be done in a configurable fashion
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
  
  get_gradient_energies(g_qm, g_mm, e);
}

void QM_QChem::get_gradient_energies(std::vector<double> &g_qm,
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

std::string QM_QChem::get_qcprog(void){
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

std::string QM_QChem::get_qcscratch(void){
  char * pwd = std::getenv("PWD");
  std::string scratch_path = std::string(pwd) + "/GQSH.sav";
  
  char * qc_str = std::getenv("QCSCRATCH");
  if (nullptr == qc_str){
    std::cerr << "Warning, $QCSCRATCH not set; using " << scratch_path << std::endl;
  }
  else if (qc_str[0] != '/'){
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

void QM_QChem::parse_mm_gradient(std::vector<double> &g_mm){
  /*
    efield.dat uses atomic units; its format is:

    Ex(mm1) Ey(mm1) Ez(mm1)
    ...
    Ex(mmN) Ey(mmN) Ez(mmN)
    Ex(qm1) Ey(qm1) Ez(qm1)
    ...
    Ex(qmN) Ey(qmN) Ez(qmN)

    where Ea(u) is the component of the electric field in the a
    direction at u and mmi and qmi are the coordinates of the ith mm
    and qm atoms respectively.
  */

  /*
    FIXME: This will return a *force* rather than a gradient unless it
    is the negative of the field, which it may be. We need to return a
    gradient.
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
}

void QM_QChem::parse_qm_gradient(std::vector<double> &g_qm,
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
}


void QM_QChem::exec_qchem(void){
  std::string cmd = "cd " + qc_scratch_directory + "; " + qc_executable + " " + qc_input_file + " " + qc_scratch_directory;
  std::cout << cmd << std::endl;
  int status = std::system(cmd.c_str());
  if (status){
    throw std::runtime_error("Q-Chem could not be called or exited abnormally.");
  }
  first_call = false;
}

void QM_QChem::write_gradient_job(void){
  std::ofstream ifile(qc_scratch_directory + "/" + qc_input_file);
  ifile.setf(std::ios_base::fixed, std::ios_base::floatfield);
  ifile.precision(std::numeric_limits<double>::digits10);

  write_molecule_section(ifile);

  /*
    NOTE: require 
      input_bohr = true
    since all inputs/outputs from this module must be in atomic units
  */

  ifile << R"(
$rem
  jobtype          force
  method           hf
  basis	           6-31+G*
  sym_ignore       true
  qm_mm            true     # external charges in NAC; generate efield.dat
  qmmm_ext_gifs    1
  input_bohr       true
)";

  
  if (! first_call ){
    ifile << "  scf_guess        read" << std::endl;
  }
  
  ifile << "$end" << std::endl;

  ifile.close();
}

void QM_QChem::write_molecule_section(std::ostream &ifile){
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
  ifile << "$molecule" << std::endl;
  ifile << qm_charge << " " << qm_multiplicity << std::endl;

  for (size_t i = 0; i < NQM; i++){
    ifile << atomids[i]      << " ";        // id
    ifile << crd_qm[i*3 + 0] << " ";        // x
    ifile << crd_qm[i*3 + 1] << " ";        // y
    ifile << crd_qm[i*3 + 2] << std::endl;  // x
  }
  ifile <<  "$end" << std::endl;

  if (NMM > 0){
    ifile << std::endl << "$external_charges" << std::endl;    
    for (size_t i = 0; i < NMM; i++){
      ifile << crd_mm[i*3 + 0] << " ";        // x
      ifile << crd_mm[i*3 + 1] << " ";        // y
      ifile << crd_mm[i*3 + 2] << " ";        // z
      ifile << chg_mm[i]       << std::endl;  // charge
    }
    ifile << "$end" << std::endl;
  }
}

int QM_QChem::readQFMan(int filenum, std::vector<double> &v){
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
