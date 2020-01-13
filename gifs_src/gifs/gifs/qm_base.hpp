#ifndef GIFS_SH_QM_BASE_H
#define GIFS_SH_QM_BASE_H

#include <vector>

QMInterface::QMInterface(size_t nqm, std::vector<int> &qmid) 
    : NQM{nqm}, NMM{0}, crd_qm{}, atomids{qmid}, crd_mm{}, chg_mm{}, first_call{true} 
{
      crd_qm.resize(3 * nqm);
};

template<typename T>
void
QMInterface::update(const T* crdqm, size_t nmm, const T* crdmm, const T* chgmm) {

    if (chg_mm.size() < nmm) {
        chg_mm.resize(nmm);
        crd_mm.resize(nmm*3);
    }    

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

};

#endif
