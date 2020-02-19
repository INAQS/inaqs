#include "linkatoms.hpp"


static 
LinkAtoms 
LinkAtoms::with_const_factors(arma::umat global_idx, arma::vec factors)
{
    return LinkAtoms(global_idx, factors, arma::vec{});
};
//
static 
LinkAtoms 
LinkAtoms::with_const_length(arma::umat global_idx, arma::vec dist)
{
    arma::vec fac(dist.size());
    return LinkAtoms(global_idx, fac, dist);
};
// assume global coords
template<typename T>
void 
LinkAtoms::update_crd(T* syscrd) {
    if (!dist.empty()) {
        update_factors(syscrd);
    }
    auto itr = crd.begin();
    for (size_t idx=0; idx<nlink; ++idx) {
        auto col = global_idx.col(idx);
        arma::uword qmid = col[0]*3;
        arma::uword mmid = col[1]*3;
        double fac = factors[idx];

        for(size_t ixyz=0; ixyz<3; ++ixyz) {
            *itr++ = get_la_crd(fac, syscrd[qmid + ixyz], syscrd[mmid + ixyz]);
        }
    }
};
// assume local coords, qm and mm section
template<typename T>
void 
LinkAtoms::update_crd(T* qmcrd, T* mmcrd) {
    if (!dist.empty()) {
        update_factors(qmcrd, mmcrd);
    }

    auto itr = crd.begin();
    for (size_t idx=0; idx<nlink; ++idx) {
        auto col = local_idx.col(idx);
        arma::uword qmid = col[0]*3;
        arma::uword mmid = col[1]*3;
        double fac = factors[idx];

        for(size_t ixyz=0; ixyz<3; ++ixyz) {
            *itr++ = get_la_crd(fac, qmcrd[qmid + ixyz], mmcrd[mmid + ixyz]);
        }
    }
};
    
template<typename T>
void 
LinkAtoms::update_factors(T* syscrd) {
    auto fac_itr = factors.begin();
    auto dist_itr = dist.begin();
    for (size_t idx=0; idx<nlink; ++idx) {
        auto col = local_idx.col(idx);
        arma::uword qmid = col[0]*3;
        arma::uword mmid = col[1]*3;
        double R = 0;
        for(size_t ixyz=0; ixyz<3; ++ixyz) {
            R += pow(syscrd[qmid + ixyz] - syscrd[mmid + ixyz], 2);
        }
            *fac_itr++ = *dist_itr++/std::sqrt(R);
    }
};
// assume local coords, qm and mm section
template<typename T>
void 
LinkAtoms::update_factors(T* qmcrd, T* mmcrd) {
    auto fac_itr = factors.begin();
    auto dist_itr = dist.begin();
    for (size_t idx=0; idx<nlink; ++idx) {
        auto col = local_idx.col(idx);
        arma::uword qmid = col[0]*3;
        arma::uword mmid = col[1]*3;
        double R = 0;
        for(size_t ixyz=0; ixyz<3; ++ixyz) {
            R += pow(qmcrd[qmid + ixyz] - mmcrd[mmid + ixyz], 2);
        }
        *fac_itr++ = *dist_itr++/std::sqrt(R);
    }
};

template void LinkAtoms::update_factors(double*);
template void LinkAtoms::update_factors(float*);
template void LinkAtoms::update_factors(double*, double*);
template void LinkAtoms::update_factors(float*, float*);
template void LinkAtoms::update_crd(double*);
template void LinkAtoms::update_crd(float*);
template void LinkAtoms::update_crd(double*, double*);
template void LinkAtoms::update_crd(float*, float*);
