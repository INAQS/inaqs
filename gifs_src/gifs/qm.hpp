#ifndef __QM_HPP
#define __QM_HPP

#include "properties.hpp"
#include <vector>

void p3vec(std::vector<double> vec);

class QMInterface
{
public:
  QMInterface(int nqm, std::vector<int> &qmid);
    void get_properties(PropMap &props);

    void update(std::vector<double> &crdqm,
		std::vector<double> &crdmm,
		std::vector<double> &chgmm);  

protected:
  void get_gradient(std::vector<double> &g_qm, std::vector<double> &g_mm);
  void write_gradient_job(std::string fname);
  void exec_qchem(std::string qcprog, std::string ifname, std::string savdir);
  void parse_qm_gradient(std::string savdir, std::vector<double> &g_qm);
    int NQM;             // const, actually NQM+NLink
    //  fixed size
    std::vector<int> atomids;        // NQM
    std::vector<double> crd_qm;      // NQM*3
    // flexible
    int NMM;
    std::vector<double> crd_mm;      // NMM*3
    std::vector<double> chg_mm;      // NMM

};

#endif //__QM_HPP
