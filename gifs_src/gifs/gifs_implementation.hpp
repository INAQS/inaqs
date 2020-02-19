#ifndef GIFS_SH_GIFS_IMPLEMENTATION_HPP
#define GIFS_SH_GIFS_IMPLEMENTATION_HPP

#include <vector>
#include "bomd.hpp"
#include "conversion.hpp"
#include "linkatoms.hpp"
#include <armadillo>

/* Virtual base class, defining all operations that 
 * should be called from the MD software */

// factory to create new bomd instance, based on an input file!


class GifsImpl
{
public: 
    // creation
    static GifsImpl* get_instance(size_t nqm, const int * qmid);
    static GifsImpl* get_instance();
    // methods
    template<typename T>
    T update_gradient(const T* in_qm_crd, const size_t* local_index, 
		      size_t in_nmm, const T* in_mm_crd, const T* in_mm_chg, 
		      T* in_qm_frc, T* in_mm_frc);
    //
    template<typename T> 
    void rescale_velocities(T* total_gradient, T* masses, T* velocities);
    //
    static inline bool instance_exists() noexcept {
        if (impl == nullptr) {
            return false;
        }
        return true;
    }
private:
    GifsImpl(size_t nqm, const int * qmids);
    // private destructor
    static void destory_instance();
    //
    static GifsImpl* impl; 
    // Local data!
    arma::uword nqm;
    // 
    arma::uvec qm_atomicnumbers{};
    // global index matrix
    arma::uvec qm_index{};
    arma::uvec mm_index{};
    // coordinates
    arma::mat qm_crd{};
    arma::mat mm_crd{};
    // parameters
    arma::vec mm_chg{};
    // forces
    arma::mat qm_frc{};
    arma::mat mm_frc{};
    // Classes
    Conversion* conv{nullptr};
    BOMD* bomd{nullptr};
    LinkAtoms* las{nullptr};
};


class Gifs
{
public:
    // create instance, can only be called once!
    explicit Gifs(size_t nqm, const int * qmid) {
        if (!GifsImpl::instance_exists()) {
            impl = GifsImpl::get_instance(nqm, qmid);   
        } else {
            impl = GifsImpl::get_instance();   
        }
    }
    // get a local handle to the interface 
    explicit Gifs() {
        impl = GifsImpl::get_instance();   
    }
    //
    template<typename T>
    inline
    T update_gradient(const T* qm_crd, const size_t* local_index, 
                      size_t nmm, const T* mm_crd, const T* mm_chg, 
                      T* qm_frc, T* mm_frc) {
        impl->update_gradient(qm_crd, local_index, 
                              nmm, mm_crd, mm_chg, 
                              qm_frc, mm_frc);
    }
    //
    template<typename T> 
    inline
    void rescale_velocities(T* total_gradient, T* masses, T* velocities) {
        impl->rescale_velocities(total_gradient, masses, velocities);
    }
    //
private:
    GifsImpl* impl{nullptr}; 
};

#endif // GIFS_SH_GIFS_CORE_H
