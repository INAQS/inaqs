#ifndef GIFS_LINKATOMS_H
#define GIFS_LINKATOMS_H

#include <armadillo>

class LinkAtoms
{
public:

    static LinkAtoms with_const_factors(const arma::umat global_idx, arma::vec factors)
    {
        return LinkAtoms(global_idx, factors, arma::vec{});
    };
    //
    static LinkAtoms with_const_length(const arma::umat global_idx, arma::vec dist)
    {
        arma::vec fac(dist.size());
        return LinkAtoms(global_idx, fac, dist);
    };
    //
    arma::mat& get_frc() { return frc; }
    arma::mat& get_crd() { return crd; }
    //
    inline void set_local_idx(const arma::uword* indices) {
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
              global_idx{in_idxs}, factors{in_factors}, 
              dist{in_dist}, nlink{in_factors.size()}, 
              crd(3, nlink), frc(3, nlink), local_idx(2, nlink)
              {}
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

#endif // GIFS_LINKATOMS_H
