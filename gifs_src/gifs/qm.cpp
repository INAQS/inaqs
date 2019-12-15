#include <vector>
#include <iostream>
#include <stdexcept>
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
  
  for(size_t i = 0; i < g_qm.size(); i++){
    g_qm[i] = 0;
  }

  for(size_t i = 0; i < g_mm.size(); i++){
    g_mm[i] = 0;
  }
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
