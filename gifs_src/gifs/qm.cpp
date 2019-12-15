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

void p3vec(std::vector<double> vec){
  if (0 != vec.size()%3){
    throw std::invalid_argument("Vector not in R3!");
  }

  double x,y,z;
  
  for (size_t i = 0; i < vec.size()/3; i++){
    x = vec[i*3 + 0];
    y = vec[i*3 + 1];
    z = vec[i*3 + 2];
    printf(FMT "\t" FMT "\t" FMT "\n", x, y, z);
  }
}

QMInterface::QMInterface(int nqm, std::vector<int> &qmid){
  NQM = nqm;
  atomids = qmid;
}

void QMInterface::update(std::vector<double> &crdqm,
			 std::vector<double> &crdmm,
			 std::vector<double> &chgmm){
  crd_qm = crdqm;
  NMM = chgmm.size();
  crd_mm = crdmm;
  chg_mm = chgmm;
}

void QMInterface::get_properties(PropMap &props){
  std::vector<double> g_qm, g_mm;

  g_qm=props.get(QMProperty::qmgradient);
  g_mm=props.get(QMProperty::mmgradient);

  get_gradient(g_qm, g_mm);
  //FIXME: the gradient doesn't make it above this call; some referenceing issue???
}

void QMInterface::get_gradient(std::vector<double> &g_qm,
			       std::vector<double> &g_mm){
  std::string ifname = "GQSH.in";
  std::string qcprog = "../../qcprog.exe";
  std::string savdir = "./GQSH.sav";
  write_gradient_job(ifname);
  exec_qchem(qcprog, ifname, savdir);
  parse_qm_gradient(savdir, g_qm);
}

void QMInterface::parse_qm_gradient(std::string savdir,
				    std::vector<double> &g_qm){
  std::string gfile = savdir + "/" + "GRAD";
  std::ifstream ifile(gfile);
  std::string line;

  bool gradient = false;
  bool energy = false;
  int i = 0;
  double E = 0;
  
  while( getline(ifile, line) && i < NQM){
    std::stringstream s(line);
    if (energy){
      s>>E;
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
}

void QMInterface::exec_qchem(std::string qcprog,
			     std::string ifname,
			     std::string savdir){
  std::string cmd = qcprog + " " + ifname + " " + savdir;
  std::system(cmd.c_str());
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
$end
)" << std::endl;
  
  ifile << "$molecule \n0 1" << std::endl;
  
  double x, y, z;
  {
    int id;
    for (int i = 0; i < NQM; i++){
      x = crd_qm[i*3 + 0];
      y = crd_qm[i*3 + 1];
      z = crd_qm[i*3 + 2];
      id = atomids[i];
      
      ifile << id << " ";
      ifile << x << " ";
      ifile << y << " ";
      ifile << z << std::endl;
    }
  }
  ifile <<  "$end" << std::endl << std::endl;

  ifile << "$external_charges" << std::endl;

  {
    double chg;
    for (int i = 0; i < NQM; i++){
      x = crd_mm[i*3 + 0];
      y = crd_mm[i*3 + 1];
      z = crd_mm[i*3 + 2];
      chg = chg_mm[i];
      
      ifile << x << " ";
      ifile << y << " ";
      ifile << z << " ";
      ifile << chg << std::endl;
    }
  }
  ifile << "$end" << std::endl;

  ifile.close();
}

int main(void){
  std::vector<double> qm_crd, mm_crd, mm_chg;
  int nqm = 4;
  int nmm = 2;

  std::vector<int> qmids = {9, 9, 9, 9};
  
  QMInterface* qm = new QMInterface(nqm, qmids);
  
  qm_crd = { 2.1194186904,    -0.1478987347,     0.0000000000,
	     1.4498754181,     1.0202212338,     0.0000000000,
	     -1.4498754181,   -1.0202212338,     0.0000000000,
	     -2.1194186904,    0.1478987347,     0.0000000000};
  
  mm_crd = {0.00000, 0.00000,  1.00000,
	    0.00000, 0.00000, -1.00000};

  mm_chg = {0.1, -0.1};
  
  std::cout << "qm_crd" << std::endl; p3vec(qm_crd);
  std::cout << "mm_crd" << std::endl; p3vec(mm_crd);

  qm->update(qm_crd, mm_crd, mm_chg);
  
  std::vector<double> g_qm, g_mm;
  g_qm.resize((size_t) 3 * nqm);
  g_mm.resize((size_t) 3 * nmm);
  
  PropMap props;
  props.emplace(QMProperty::qmgradient, g_qm);
  props.emplace(QMProperty::mmgradient, g_mm);
  qm->get_properties(props);
  
  std::cout << "g_qm" << std::endl; p3vec(g_qm);
  std::cout << "g_mm" << std::endl; p3vec(g_mm);

  return 0;
}
