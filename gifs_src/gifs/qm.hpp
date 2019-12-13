#ifndef __LAUNCH_HPP
#define __LAUNCH_HPP
#include <unordered_map>
#include <vector>

enum QMProps{
	     qmgradient,
	     mmgradient,
	     energy,
};

using PropMap = std::unordered_map<QMProps, std::vector<double>&>;

class QM{
public:
  QM(int nqm, int nmm,
     std::vector<int> qmid,
     std::vector<double> mmcharge);
  
  void get_properties(PropMap &props,
		      std::vector<double> &qm_crd,
		      std::vector<double> &mm_crd);  
protected:
  int n_qm;
  int n_mm;
  
  std::vector<int> qm_id;
  std::vector<double> mm_charge;
};

#endif
