#include <vector>
#include <iostream>
#include <stdexcept>
#include "launch.hpp"

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

std::vector<double>& get(PropMap m, QMProps p){
  auto itr = m.find(p);
  if(itr == m.end()){
    throw std::invalid_argument("Bad Key!");
  }
  else return itr->second;
}

QM::QM(int nqm, int nmm,
       std::vector<int> qmid,
       std::vector<double> mmcharge){
  n_qm = nqm;
  n_mm = nmm;
  
  qm_id = qmid;
  mm_charge = mmcharge;
}
  
void QM::get_properties(PropMap &props,
			std::vector<double> &qm_crd,
			std::vector<double> &mm_crd){
  std::vector<double> g_qm, g_mm;

  g_qm=get(props, qmgradient);
  g_mm=get(props, mmgradient);
  
  for(size_t i = 0; i < g_qm.size(); i++){
    g_qm[i] = 0;
  }

  for(size_t i = 0; i < g_mm.size(); i++){
    g_mm[i] = 0;
  }
}

int main(void){
  std::vector<double> qm_crd, mm_crd;
  int nqm = 4;
  int nmm = 2;
  
  QM* qm = new QM(nqm, nmm, std::vector<int> {9, 9, 9, 9},
  		  std::vector<double> {0.1, -0.1});

  qm_crd = { 2.1194186904,    -0.1478987347,     0.0000000000,
	     1.4498754181,     1.0202212338,     0.0000000000,
	     -1.4498754181,   -1.0202212338,     0.0000000000,
	     -2.1194186904,    0.1478987347,     0.0000000000};
  
  mm_crd = {0.00000, 0.00000,  1.00000,
	    0.00000, 0.00000, -1.00000};
  
  std::cout << "qm_crd" << std::endl; p3vec(qm_crd);
  std::cout << "mm_crd" << std::endl; p3vec(mm_crd);

  std::vector<double> g_qm, g_mm;
  g_qm.resize((size_t) 3 * nqm);
  g_mm.resize((size_t) 3 * nmm);

  PropMap props = {{qmgradient, g_qm}, {mmgradient, g_mm}};

  qm->get_properties(props, qm_crd, mm_crd);
  
  std::cout << "g_qm" << std::endl; p3vec(g_qm);
  std::cout << "g_mm" << std::endl; p3vec(g_mm);
  
  return 0;
}
