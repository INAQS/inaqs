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
