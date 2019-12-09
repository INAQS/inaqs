#ifndef GIFS_SH_GIFS_CORE_H
#define GIFS_SH_GIFS_CORE_H

#include <string>

namespace gifs 
{
/* Virtual base class, defining all operations that 
 * should be called from the MD software */

// factory to create new bomd instance, based on an input file!
BOMD* get_bomd_instance(std::string fname) {
    // store objects in a global map, so that
    // the addition of user defined BOMD codes is possible
};


class Gifs
{
    // create instance, can only be called once!
    explicit Gifs(std::string fname) {
        impl = GifsImpl::get_instance(fname);   
    }
    // get a local handle to the interface 
    explicit Gifs() {
        impl = GifsImpl::get_instance();   
    }

    template<typename T>
    inline
    T get_gradient(T* qm_crd, T* mm_crd, T* mm_chg, T* qm_gradient, T* mm_gradient)
    {
        impl->get_gradient(qm_crd, mm_crd, mm_chg, qm_gradient, mm_gradient);
    }

    template<typename T>
    inline 
    T rescale_velocities(T* total_gradient, T* masses, T* velocities)
    {
        impl->rescale_velocities(total_gradient, masses, velocities);
    }
private:
    GifsImpl* impl; 
};


class GifsImpl
{
public: 

    template<typename T>
    inline
    T get_gradient(T* qm_crd, T* mm_crd, T* mm_chg, T* qm_gradient, T* mm_gradient)
    {
        bomd->get_gradient(qm_crd, mm_crd, mm_chg, qm_gradient, mm_gradient);
    };

    template<typename T>
    inline 
    T rescale_velocities(T* total_gradient, T* masses, T* velocities)
    {
        bomd->rescale_velocities(total_gradient, masses, velocities);
    };

    // creation
    static 
    GifsImpl* get_instance(std::string fname)
    {
        // actually create instance, and register 
        // its destructor atexit!
        if (impl != nullptr)
            throw std::Exception("GIFS object was already created");
        impl = new GifsImpl{fname};
        //
        atexit(GifsImpl::destory_instance);
        return impl;
    };

    static 
    GifsImpl* get_instance() {
        if (impl == nullptr) {
            throw std::Exception("GIFS object was not created, yet");
        }
        return impl;
    };

private:
    //
    GifsImpl(std::string fname) : bomd{get_bomd(fname)} {};
    //
    static destory_instance() {
        if (impl != nullptr) {
            delete impl->bomd;
            delete impl;
        }
    }
    //
    static GifsImpl* impl; 
    static int icount; 
    //
    BOMD* bomd{nullptr};
};

GifsImpl::impl = nullptr;
GifsImpl::icount = 0;

} // end namespace gifs

#endif // GIFS_SH_GIFS_CORE_H
