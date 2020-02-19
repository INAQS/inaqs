#include "gifs_implementation.hpp"
#include <stdlib.h>
#include <vector>
#include <algorithm>

// creation
GifsImpl* GifsImpl::get_instance(size_t nqm, const int * qmid)
{
// actually create instance, and register its destructor atexit!
    if (instance_exists()) {
        // throw std::Exception("GIFS object was already created");
    }
    //
    impl = new GifsImpl(nqm, qmid);
    //
    atexit(destory_instance);
    return impl;
};
//
GifsImpl* GifsImpl::get_instance() {
    if (!instance_exists()) {
        //   throw std::Exception("GIFS object was not initilized, yet!");
    }
    return impl;
};
//
template<typename T>
T GifsImpl::update_gradient(const T* in_qm_crd, const size_t* local_index, 
                            size_t in_nmm, const T* in_mm_crd, const T* in_mm_chg, 
                            //output
                            T* in_qm_frc, T* in_mm_frc)
{
    //
    las->set_local_idx(local_index);
    // set mm size
    arma::resize(mm_crd, 3, nmm);
    arma::resize(mm_frc, 3, nmm);
    arma::resize(mm_chg, nmm);
    // update charge
    for (arma::uword idx=0; idx<nmm; ++idx) {
        mm_chg[idx] = in_mm_chg[idx];
    };
    // update linkatom coordinates and set factors etc.
    las->update_crd(in_qm_crd, in_mm_crd);
    // copy coordinates 
    size_t nqm_withoutlink = nqm-las->nlink;
    conv->transform_coords_to_au(in_qm_crd, in_qm_crd+3*nqm_withoutlink, qm_crd.begin());
    conv->transform_coords_to_au(in_mm_crd, in_mm_crd+3*nmm, mm_crd.begin()); 
    // linkatoms coords
    const auto& la_crd = las->get_crd();
    conv->transform_coords_to_au(la_crd.begin(), la_crd.end(), qm_crd.begin()+nqm_withoutlink*3);
    // Compute Forces etc.
    double energy = bomd->update_gradient();
    // update linkatom
    auto& la_frc = las->get_frc();
    conv->transform_coords_to_au(qm_frc.begin()+3*nqm_withoutlink, qm_frc.end(), la_frc.begin());
    // transform back forces
    conv->transform_coords_to_au(qm_frc.begin(), qm_frc.end()-3*las->nlink, in_qm_frc);
    conv->transform_coords_to_au(mm_frc.begin(), mm_frc.end(), in_mm_frc);
    // copy linkatom forces
    std::copy(la_frc.begin(), la_frc.end(), in_qm_frc+3*nqm_withoutlink);
    // return energy
    return conv->energy_from_au(energy);
};
//
template<typename T> 
void GifsImpl::rescale_velocities(T* total_gradient, T* masses, T* velocities)
{
    bomd->rescale_velocities(total_gradient, masses, velocities);
};
//
GifsImpl::GifsImpl(size_t in_nqm, const int * ian)
    : nqm{in_nqm}, qm_atomicnumbers(in_nqm), qm_index(in_nqm),
      mm_index(0), qm_crd(3, in_nqm), mm_crd(3, 1), mm_chg(1),
      qm_frc(3, in_nqm), mm_frc(3, 1)
{
    // copy atomic numbers
    for (arma::uword idx=0; idx<in_nqm; ++idx) {
        qm_atomicnumbers[idx] = ian[idx];
    }
    // should be setup by initialization
    const double mass_unit = 1.660538921e-27; 
    const double length_unit = 1e-9; 
    const double time_unit = 1e-12;
    //
    //bomd = new BOMD(in_nqm, in_qmid);
    BOMD = new BOMD(qm_atomicnumbers,
		    qm_crd, 
		    mm_crd, 
		    mm_chg, 
		    qm_grd,  //FIXME: MFSJM: Where are these gradients living?
		    mm_grd);
    // 
    conv = Conversion::from_elementary(mass_unit, length_unit, time_unit);
    // no link atoms
    arma::umat global_idx(0, 0);
    arma::vec factors(0);
    las = LinkAtoms::with_const_factors(global_idx, factors);
};
//
void GifsImpl::destory_instance() {
    if (instance_exists()) {
        delete impl->bomd;
        delete impl->conv;
        delete impl->las;
        delete impl;
    }
}
//
GifsImpl* GifsImpl::impl = nullptr;
//
template void GifsImpl::rescale_velocities(double* total_gradient, double* masses, double* velocities);
template void GifsImpl::rescale_velocities(float* total_gradient, float* masses, float* velocities);
template double GifsImpl::get_gradient(const double* qm_crd, size_t nmm, const double* mm_crd, const double* mm_chg, double* f, double* fshift);
template float GifsImpl::get_gradient(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift);
