#ifndef GIFS_SH_GIFS_CORE_H
#define GIFS_SH_GIFS_CORE_H

#include <vector>
#include "bomd.hpp"

/* Virtual base class, defining all operations that 
 * should be called from the MD software */

// factory to create new bomd instance, based on an input file!


class GifsImpl
{
public: 

    template<typename T> inline
    T get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift);

    template<typename T> inline
    void rescale_velocities(T* total_gradient, T* masses, T* velocities);

    // creation
    static GifsImpl* get_instance(size_t nqm, std::vector<int>& qmid);

    static GifsImpl* get_instance();
    
    static inline bool instance_exists() noexcept;

private:
    //
    GifsImpl(size_t nqm, std::vector<int>& qmid);
    //
    static void destory_instance();
    //
    static GifsImpl* impl; 
    //
    BOMD* bomd{nullptr};
};


class Gifs
{
public:
    // create instance, can only be called once!
    explicit Gifs(size_t nqm, std::vector<int>& qmid);
    // get a local handle to the interface 
    explicit Gifs();

    template<typename T> inline
    T get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift);

    template<typename T> inline
    void rescale_velocities(T* total_gradient, T* masses, T* velocities);

private:
    GifsImpl* impl{nullptr}; 
};

#endif // GIFS_SH_GIFS_CORE_H
