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

QMInterface::QMInterface(size_t nqm, std::vector<int> &qmid){
  NQM = nqm;
  atomids = qmid;
  crd_qm.reserve(nqm);
}

void QMInterface::update(std::vector<double> &crdqm,
			 std::vector<double> &crdmm,
			 std::vector<double> &chgmm){
  crd_qm = crdqm;
  NMM = chgmm.size();
  crd_mm = crdmm;
  chg_mm = chgmm;
}

void QMInterface::update(const float* crdqm, size_t nmm, const float* crdmm, const float* chgmm)
{
  NMM = nmm;
  if (nmm < crd_mm.capacity()) {
        crd_mm.resize(nmm);
        chg_mm.resize(nmm);
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

  write_gradient_job(ifname);
  exec_qchem(qcprog, ifname, savdir);
  parse_qm_gradient(savdir, g_qm, e);
}

std::string QMInterface::get_qcprog(void){
  char * qc_str = std::getenv("QC");
  std::string default_path = "../../exe/qcprog.exe";
  if (nullptr == qc_str){
    std::cerr << "Warning, $QC not set; using " << default_path << std::endl;
    return default_path;
  }
  else{
    return std::string(qc_str) + "/exe/qcprog.exe";
  }
}

void QMInterface::parse_qm_gradient(std::string savdir,
				    std::vector<double> &g_qm,
				    std::vector<double> &e){
  std::string gfile = savdir + "/" + "GRAD";
  std::ifstream ifile(gfile);
  std::string line;

  bool gradient = false;
  bool energy = false;
  size_t i = 0;
  
  while( getline(ifile, line) && i < NQM){
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

  /*
    Unit conversion: Q-Chem outputs (??? need to confirm) in
    Hartree/Angstrom; our interface expects Hartree/Bohr
  */
  {
    const double ang2bohr = 0.529177249;
    for (auto& v : g_qm){
      v *= ang2bohr;
    }
  }
}

void QMInterface::exec_qchem(std::string qcprog,
			     std::string ifname,
			     std::string savdir){
  std::string cmd = qcprog + " " + ifname + " " + savdir;
  int status = std::system(cmd.c_str());
  if (status){
    throw std::runtime_error("Q-Chem could not be called or exited abnormally.");
  }
  first_call = false;
}

void QMInterface::write_gradient_job(std::string fname){
  std::ofstream ifile(fname);

  ifile << R"(
$rem
  jobtype          force
  method           hf
  basis	           6-31+G*
  sym_ignore       true
  qm_mm            true     # external charges in NAC; generate efield.dat
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

  
  
  ifile << "$molecule \n0 1" << std::endl;
  /*
    format of $molecule section:
    [atomic-number] [x-coord] [y-coord] [x-coord]
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
    
    for (size_t i = 0; i < NQM; i++){
      ifile << crd_mm[i*3 + 0] << " ";        // x
      ifile << crd_mm[i*3 + 1] << " ";        // y
      ifile << crd_mm[i*3 + 2] << " ";        // z
      ifile << chg_mm[i]       << std::endl;  // charge
    }
    ifile << "$end" << std::endl;
  }

  ifile.close();
}

