#ifndef GIFS_SH_QM_BASE_H
#define GIFS_SH_QM_BASE_H

#include <vector>

template<typename T>
void
QMInterface::update(const T* crdqm, size_t nmm, const T* crdmm, const T* chgmm) {
  NMM = nmm;
  if (chg_mm.size() < nmm){
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

#endif
