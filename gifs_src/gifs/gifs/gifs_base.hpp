#ifndef GIFS_SH_GIFS_BASE_H
#define GIFS_SH_GIFS_BASE_H

#include <stdlib.h>

/* Virtual base class, defining all operations that 
 * should be called from the MD software */

// factory to create new bomd instance, based on an input file!


template<typename T>
inline
T GifsImpl::get_gradient(T* qm_crd, T* mm_crd, T* mm_chg, T* qm_gradient, T* mm_gradient)
{
    bomd->get_gradient(qm_crd, mm_crd, mm_chg, qm_gradient, mm_gradient);
};

template<typename T>
inline 
T GifsImpl::rescale_velocities(T* total_gradient, T* masses, T* velocities)
{
    bomd->rescale_velocities(total_gradient, masses, velocities);
};

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
}

template<typename T>
inline
T Gifs::get_gradient(T* qm_crd, T* mm_crd, T* mm_chg, T* qm_gradient, T* mm_gradient)
{
    impl->get_gradient(qm_crd, mm_crd, mm_chg, qm_gradient, mm_gradient);
}

template<typename T>
inline 
T Gifs::rescale_velocities(T* total_gradient, T* masses, T* velocities)
{
    impl->rescale_velocities(total_gradient, masses, velocities);
}


#endif // GIFS_SH_GIFS_CORE_H
