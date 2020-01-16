#include "gifs_implementation.hpp"
#include <stdlib.h>
#include <vector>

/* Virtual base class, defining all operations that 
 * should be called from the MD software */

// factory to create new bomd instance, based on an input file!
template<typename T>
T GifsImpl::get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift)
{
    return bomd->get_gradient(qm_crd, nmm, mm_crd, mm_chg, f, fshift);
};

template double GifsImpl::get_gradient(const double* qm_crd, size_t nmm, const double* mm_crd, const double* mm_chg, double* f, double* fshift);
template float GifsImpl::get_gradient(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift);

template<typename T> 
void GifsImpl::rescale_velocities(T* total_gradient, T* masses, T* velocities)
{
    bomd->rescale_velocities(total_gradient, masses, velocities);
};

template void GifsImpl::rescale_velocities(double* total_gradient, double* masses, double* velocities);
template void GifsImpl::rescale_velocities(float* total_gradient, float* masses, float* velocities);

    // creation
GifsImpl* GifsImpl::get_instance(size_t nqm, std::vector<int>& qmid)
{
// actually create instance, and register 
// its destructor atexit!
/*
    if (impl != nullptr)
    // throw std::Exception("GIFS object was already created");
*/
    impl = new GifsImpl(nqm, qmid);
    //
    atexit(destory_instance);
    return impl;
};

GifsImpl* GifsImpl::get_instance() {
    if (!instance_exists()) {
        //   throw std::Exception("GIFS object was not initilized, yet!");
    }
    return impl;
};
    

inline bool GifsImpl::instance_exists() noexcept {
    if (impl == nullptr) {
        return false;
    }
    return true;
}


GifsImpl::GifsImpl(size_t nqm, std::vector<int>& qmid) : bomd{new BOMD(nqm, qmid)} {};
//
void GifsImpl::destory_instance() {
    if (impl != nullptr) {
        delete impl->bomd;
        delete impl;
    }
}


GifsImpl* GifsImpl::impl = nullptr;


// create instance, can only be called once!
Gifs::Gifs(size_t nqm, std::vector<int>& qmid) {
    if (!GifsImpl::instance_exists()) {
        impl = GifsImpl::get_instance(nqm, qmid);   
    } else {
        impl = GifsImpl::get_instance();   
    }
}


// get a local handle to the interface 
Gifs::Gifs() {
    impl = GifsImpl::get_instance();   
};

template<typename T>
T Gifs::get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift)
{
    return impl->get_gradient(qm_crd, nmm, mm_crd, mm_chg, f, fshift);
};

template double Gifs::get_gradient(const double* qm_crd, size_t nmm, const double* mm_crd, const double* mm_chg, double* f, double* fshift);
template float Gifs::get_gradient(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift);

template<typename T>
void Gifs::rescale_velocities(T* total_gradient, T* masses, T* velocities)
{
    impl->rescale_velocities(total_gradient, masses, velocities);
};

template void Gifs::rescale_velocities(double* total_gradient, double* masses, double* velocities);
template void Gifs::rescale_velocities(float* total_gradient, float* masses, float* velocities);
