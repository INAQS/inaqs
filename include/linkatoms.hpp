#ifndef GIFS_LINKATOMS_H
#define GIFS_LINKATOMS_H

#include <armadillo>

class LinkAtoms
{
public:

    // you own it
    static LinkAtoms* with_const_factors(const arma::umat global_idx, arma::vec factors);
    //
    static LinkAtoms* with_const_length(const arma::umat global_idx, arma::vec dist);
    //
    arma::mat& get_frc() { return frc; }
    arma::mat& get_crd() { return crd; }
    //
    inline void set_local_idx(const size_t* indices) {
        auto itr = local_idx.begin();
        for (size_t idx=0; idx<nlink; ++idx) {
            *itr++ = indices[idx*2];
            *itr++ = indices[idx*2+1];
        }
    };
    // assume global coords
    template<typename T>
    void update_crd(T* syscrd);
    // assume local coords!
    template<typename T>
    void update_crd(T* qmcrd, T* mmcrd);
    //
public:
    arma::uword nlink{0};
private:
    LinkAtoms(arma::umat in_idxs, arma::vec in_factors, arma::vec in_dist):
      nlink{in_factors.size()}, factors{in_factors},  dist{in_dist},
      crd(3, nlink), frc(3, nlink), global_idx{in_idxs}, local_idx(2, nlink) {}
    //
    template<typename T>
    inline double 
    get_la_crd(const double& x, const T& qmcrd, const T& mmcrd) {
        return (1.0-x) * qmcrd + x * mmcrd;
    };
    //
    template<typename T>
    void update_factors(T* qmcrd, T* mmcrd);
    //
    template<typename T>
    void update_factors(T* syscrd); 
    //
protected:
     //
    arma::vec factors{};
    arma::vec dist{};
    // 
    arma::mat crd{};
    arma::mat frc{};
    // global idx
    arma::umat global_idx{};
    // local index
    arma::umat local_idx{};
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

#endif // GIFS_LINKATOMS_H
