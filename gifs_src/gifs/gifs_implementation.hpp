#ifndef GIFS_SH_GIFS_IMPLEMENTATION_HPP
#define GIFS_SH_GIFS_IMPLEMENTATION_HPP

#include <vector>
#include "bomd.hpp"
#include "conversion.hpp"
#include "linkatoms.hpp"

/* Virtual base class, defining all operations that 
 * should be called from the MD software */

// factory to create new bomd instance, based on an input file!


class GifsImpl
{
public: 

    template<typename T> 
    inline
    T get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f);

    template<typename T> 
    inline
    void 
    rescale_velocities(T* total_gradient, T* masses, T* velocities);

    // creation
    static GifsImpl* get_instance(size_t nqm, const int * qmid);

    static GifsImpl* get_instance();
    
    static inline bool instance_exists() noexcept;


private:
    void update_coords(size_t nmm, double* crd);
    void update_coords(const double* qm_crd, const size_t nmm, const double* mm_crd);
    //
    void set_forces(double* forces);
    // Helper
    template<typename itr_qm, typename itr_mm>
    void _update_crd(size_t nmm, itr_qm qm_start, itr_qm qm_end, itr_mm mm_start, itr_mm mm_end); 
    template<typename itr_qm>//, typename itr_mm, typename itr_la>
    void _set_forces(itr_qm qm_start);
    //
    //
    GifsImpl(size_t nqm, const int * qmids);
    //
    static void destory_instance();
    //
    static GifsImpl* impl; 
    //
    BOMD* bomd{nullptr};
    // data
    // index matrix
    std::vector<int> qm_indx{};
    std::vector<int> mm_indx{};
    // coordinates
    std::vector<double> qm_crd{};
    std::vector<double> mm_crd{};
    // parameters
    std::vector<double> mm_chg{};
    // forces
    std::vector<double> qm_frc{};
    std::vector<double> mm_frc{};
    //
    SystemConversion conv{};
    // 
    LinkAtoms las{};
};


class Gifs
{
public:
    // create instance, can only be called once!
    explicit Gifs(size_t nqm, const int * qmid);
    // get a local handle to the interface 
    explicit Gifs();

    template<typename T>
    T get_gradient(const T* qm_crd, size_t nmm, const T* mm_crd, const T* mm_chg, T* f, T* fshift);

    template<typename T> inline
    void rescale_velocities(T* total_gradient, T* masses, T* velocities);

private:
    GifsImpl* impl{nullptr}; 
};

#endif // GIFS_SH_GIFS_CORE_H
