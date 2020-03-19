#include <stdlib.h>
#include <vector>
#include <algorithm>
//
#include "gifs_implementation.hpp"

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
    nmm = in_nmm;
    //
    las->set_local_idx(local_index);
    // set mm size
    mm_crd.resize(3, nmm);
    mm_frc.resize(3, nmm);
    mm_chg.resize(nmm);
    mm_index.resize(nmm);
    // update charge
    for (arma::uword idx=0; idx<nmm; ++idx) {
        mm_chg[idx] = in_mm_chg[idx];
    };
    // update linkatom coordinates and set factors etc.
    las->update_crd(in_qm_crd, in_mm_crd);
    // copy coordinates 
    size_t nqm_withoutlink = nqm-las->nlink;
    // assumes linkatoms are below coordinate section
    conv->transform_coords_to_au(in_qm_crd, in_qm_crd+3*nqm_withoutlink, qm_crd.begin());
    conv->transform_coords_to_au(in_mm_crd, in_mm_crd+3*nmm, mm_crd.begin()); 
    // linkatoms coords
    const auto& la_crd = las->get_crd();
    conv->transform_coords_to_au(la_crd.begin(), la_crd.end(), qm_crd.begin()+nqm_withoutlink*3);
    // Compute Forces etc.
    double energy = bomd->update_gradient();
    // update linkatom
    auto& la_frc = las->get_frc();
    conv->transform_au_to_forces(qm_frc.begin()+3*nqm_withoutlink, qm_frc.end(), la_frc.begin());
    // transform back forces
    conv->transform_au_to_forces(qm_frc.begin(), qm_frc.end()-3*las->nlink, in_qm_frc);
    conv->transform_au_to_forces(mm_frc.begin(), mm_frc.end(), in_mm_frc);
    // copy linkatom forces
    std::copy(la_frc.begin(), la_frc.end(), in_qm_frc+3*nqm_withoutlink);
    // return energy
    return conv->energy_from_au(energy);
};
//
template<typename T> 
void GifsImpl::rescale_velocities(T* in_grad, T* in_masses, T* in_veloc)
{
    const arma::uword ntot = nqm + nmm;
    //
    total_gradient.resize(3, ntot);
    masses.resize(ntot);
    veloc.resize(3, ntot);

    for (arma::uword i=0; i<nqm; ++i) {
        const arma::uword idx = qm_index[i];
        total_gradient(0, i) = in_grad[idx*3];
        total_gradient(1, i) = in_grad[idx*3+1];
        total_gradient(2, i) = in_grad[idx*3+2];
    }

    for (arma::uword i=0; i<nmm; ++i) {
        const arma::uword idx = mm_index[i];
        total_gradient(0, i+nqm) = in_grad[idx*3];
        total_gradient(1, i+nqm) = in_grad[idx*3+1];
        total_gradient(2, i+nqm) = in_grad[idx*3+2];
    }

    for (arma::uword i=0; i<nqm; ++i) {
        const arma::uword idx = qm_index[i];
        veloc(0, i) = in_veloc[idx*3];
        veloc(1, i) = in_veloc[idx*3+1];
        veloc(2, i) = in_veloc[idx*3+2];
    }

    for (arma::uword i=0; i<nmm; ++i) {
        const arma::uword idx = mm_index[i];
        veloc(0, i+nqm) = in_veloc[idx*3];
        veloc(1, i+nqm) = in_veloc[idx*3+1];
        veloc(2, i+nqm) = in_veloc[idx*3+2];
    }

    for (arma::uword i=0; i<nqm; ++i) {
        const arma::uword idx = qm_index[i];
        masses[i] = in_masses[idx];
    }

    for (arma::uword i=0; i<nmm; ++i) {
        const arma::uword idx = mm_index[i];
        masses[i] = in_masses[idx];
    }
    // do unit transformation
    conv->transform_veloc_to_au(veloc.begin(), veloc.end(), veloc.begin());
    conv->transform_forces_to_au(total_gradient.begin(), total_gradient.end(), total_gradient.begin());
    double e_drift = 0;
    bomd->rescale_velocities(veloc, masses, total_gradient, e_drift);
};
//
GifsImpl::GifsImpl(size_t in_nqm, const int * ian)
    : nqm{in_nqm}, nmm{0}, qm_atomicnumbers(in_nqm), qm_index(in_nqm)
{
    // set to correct size:
    qm_crd.resize(3, nqm);
    qm_frc.resize(3, nqm);
    mm_frc.resize(3, 0);
    mm_crd.resize(3, 0);
    mm_chg.resize(0);
    mm_index.resize(0);
    qm_index.resize(nqm);
    veloc.resize(3, 0);
    masses.resize(0);
    total_gradient.resize(3, 0);
    // copy atomic numbers
    for (arma::uword idx=0; idx<in_nqm; ++idx) {
        qm_atomicnumbers[idx] = ian[idx];
    }
    // should be setup by initialization
    const double mass_unit = 1.660538921e-27; 
    const double length_unit = 1e-9; 
    const double time_unit = 1e-12;
    // 
    conv = Conversion::from_elementary(mass_unit, length_unit, time_unit);
    //bomd = new BOMD(in_nqm, in_qmid);
    bomd = new BOMD(qm_atomicnumbers,
            qm_crd, 
            mm_crd, 
            mm_chg, 
            qm_frc, 
            mm_frc);
    // no link atoms
    arma::umat global_idx(0, 0);
    arma::vec factors(1);
    factors.resize(0);
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

void Gifs::update_global_index(int* indexQM, int* indexMM) {
    impl->update_global_index(indexQM, indexMM);
}

void 
GifsImpl::update_global_index(int* indexQM, int* indexMM)
{
    for (arma::uword i=0; i<nqm; ++i) {
        qm_index[i] = indexQM[i];
    };
    for (arma::uword i=0; i<nmm; ++i) {
        mm_index[i] = indexMM[i];
    };
}
//
GifsImpl* GifsImpl::impl = nullptr;
//
template void GifsImpl::rescale_velocities(double* total_gradient, double* masses, double* velocities);
template void GifsImpl::rescale_velocities(float* total_gradient, float* masses, float* velocities);
template float GifsImpl::update_gradient(const float* in_qm_crd, const size_t* local_index, size_t in_nmm, const float* in_mm_crd, const float* in_mm_chg, float* in_qm_frc, float* in_mm_frc);
template double GifsImpl::update_gradient(const double* in_qm_crd, const size_t* local_index, size_t in_nmm, const double* in_mm_crd, const double* in_mm_chg, double* in_qm_frc, double* in_mm_frc);
