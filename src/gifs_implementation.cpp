#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
//
#include "gifs_implementation.hpp"
#include "fssh.hpp"
#include "bomd_rescale.hpp"
#include "bomd_electronic.hpp"
#include "diabatic_seam_dynamics.hpp"
#include "ehrenfest.hpp"
#include "configreader.hpp"
//
// creation
GifsImpl* GifsImpl::get_instance(const char* file, size_t nqm, const int * qmid, const double mass, const double length, const double time)
{
// actually create instance, and register its destructor atexit!
    if (instance_exists()) {
        // throw std::Exception("GIFS object was already created");
    }
    //
    impl = new GifsImpl(file, nqm, qmid, mass, length, time);
    //
    atexit(destroy_instance);
    return impl;
};
//
GifsImpl* GifsImpl::get_instance() {
    if (!instance_exists()) {
      throw std::logic_error("GIFS object was not initialized, yet!");
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
    conv.transform_crd_md2au(in_qm_crd, in_qm_crd+3*nqm_withoutlink, qm_crd.begin());
    conv.transform_crd_md2au(in_mm_crd, in_mm_crd+3*nmm, mm_crd.begin());
    
    // linkatoms coords
    const auto& la_crd = las->get_crd();
    conv.transform_crd_md2au(la_crd.begin(), la_crd.end(), qm_crd.begin()+nqm_withoutlink*3);
  
    // Compute Forces etc.
    double energy = bomd->update_gradient();
    // update linkatom
    auto& la_frc = las->get_frc();

    conv.transform_gradient_au2md(qm_frc.begin()+3*nqm_withoutlink, qm_frc.end(), la_frc.begin());
    // transform back forces
    conv.transform_gradient_au2md(qm_frc.begin(), qm_frc.end()-3*las->nlink, in_qm_frc);
    conv.transform_gradient_au2md(mm_frc.begin(), mm_frc.end(), in_mm_frc);
    // copy linkatom forces
    std::copy(la_frc.begin(), la_frc.end(), in_qm_frc+3*nqm_withoutlink);
    // return energy
    return conv.energy_au2md(energy);
};

//FIXME: pretty sure this is wrong wrt link atoms and .begin() / .end() ! 
void GifsImpl::update_coords(const arma::mat & R){
  if (R.n_cols != nmm + nqm){
    throw std::logic_error("Invalid call to gifs_update_coords(); new size does not match!");
  }
  size_t nqm_withoutlink = nqm-las->nlink;
  
  // update linkatom coordinates and set factors etc.
  las->update_crd(R.begin(), R.begin() + 3*nqm_withoutlink);
  
  // assumes linkatoms are below coordinate section
  conv.transform_crd_md2au(R.begin(), R.begin()+3*nqm_withoutlink, qm_crd.begin());
  conv.transform_crd_md2au(R.begin()+3*nqm_withoutlink, R.begin() + 3*(nqm_withoutlink + nmm), mm_crd.begin());
  
  // linkatoms coords
  const auto& la_crd = las->get_crd();
  conv.transform_crd_md2au(la_crd.begin(), la_crd.end(), qm_crd.begin()+nqm_withoutlink*3);  
}

template<typename T1, typename T2>
void from_global(const T1* global, T2& local,
                 const arma::uword nqm, const arma::uvec& qm_index,
                 const arma::uword nmm, const arma::uvec& mm_index)
{
    for (arma::uword i=0; i<nqm; ++i) {
        const arma::uword idx = qm_index[i];
        local(0, i) = global[idx*3];
        local(1, i) = global[idx*3+1];
        local(2, i) = global[idx*3+2];
    }

    for (arma::uword i=0; i<nmm; ++i) {
        const arma::uword idx = mm_index[i];
        local(0, i+nqm) = global[idx*3];
        local(1, i+nqm) = global[idx*3+1];
        local(2, i+nqm) = global[idx*3+2];
    }
};

template<typename T1, typename T2>
void to_global(T1* global, const T2& local,
               const arma::uword nqm, const arma::uvec& qm_index,
               const arma::uword nmm, const arma::uvec& mm_index)
{
    for (arma::uword i=0; i<nqm; ++i) {
        const arma::uword idx = qm_index[i];
        global[idx*3]   = local(0, i);
        global[idx*3+1] = local(1, i);
        global[idx*3+2] = local(2, i);
    }

    for (arma::uword i=0; i<nmm; ++i) {
        const arma::uword idx = mm_index[i];
        global[idx*3]   = local(0, i+nqm);
        global[idx*3+1] = local(1, i+nqm);
        global[idx*3+2] = local(2, i+nqm);
    }
};

//
template<typename T> 
void GifsImpl::rescale_velocities(T total_energy, T* in_grad, T* in_masses, T* in_veloc)
{
    const arma::uword ntot = nqm + nmm;
    //
    total_gradient.resize(3, ntot);
    masses.resize(ntot);
    veloc.resize(3, ntot);
    //
    from_global(in_grad, total_gradient, nqm, qm_index, nmm, mm_index);
    // FIXME: change names for the call tree up to this point to reflect a force
    total_gradient *= -1.0; // the in_grad parameter is actually a force
    from_global(in_veloc, veloc, nqm, qm_index, nmm, mm_index);
    //
    // FIXME: rewrite FSSH to take inv masses so we don't do this conversion
    for (arma::uword i=0; i<nqm; ++i) {
        const arma::uword idx = qm_index[i];
        masses[i] = conv.mass_md2au(1.0/in_masses[idx]);
    }

    for (arma::uword i=0; i<nmm; ++i) {
        const arma::uword idx = mm_index[i];
        masses[i+nqm] = conv.mass_md2au(1.0/in_masses[idx]);
    }
    // do unit transformation
    conv.transform_veloc_md2au(veloc.begin(), veloc.end(), veloc.begin());
    conv.transform_gradient_md2au(total_gradient.begin(), total_gradient.end(), total_gradient.begin());
    total_energy = conv.energy_md2au(total_energy);
    if (bomd->rescale_velocities(veloc, masses, total_gradient, total_energy)) {
        //
        conv.transform_veloc_au2md(veloc.begin(), veloc.end(), veloc.begin());
        conv.transform_gradient_au2md(total_gradient.begin(), total_gradient.end(), total_gradient.begin());
        total_gradient *= -1.0; // transform back to a force
        //
        to_global(in_grad, total_gradient, nqm, qm_index, nmm, mm_index);
        to_global(in_veloc, veloc, nqm, qm_index, nmm, mm_index);
    }
};
//
std::unique_ptr<BOMD>
GifsImpl::select_bomd(ConfigBlockReader& reader, FileHandle& fh,
             arma::uvec& atomicnumbers,
             arma::mat& qm_crd, 
             arma::mat& mm_crd, 
             arma::vec& mm_chg, 
             arma::mat& qm_grd,
             arma::mat& mm_grd) 
{
    std::unique_ptr<BOMD> bomd;
    std::string runtype;
    reader.get_data("runtype", runtype);
    if (runtype == "bomd") {
      bomd.reset(new BOMD(qm_grd, mm_grd));
    }
    else if (runtype == "fssh") {
      bomd.reset(new FSSH(qm_grd, mm_grd));
    }
    else if (runtype == "bomd-rescale") {
      bomd.reset(new RescaleBomd(qm_grd, mm_grd));
    }
    else if (runtype == "bomd-print") {
      bomd.reset(new PrintBomd(qm_grd, mm_grd));
    }
    else if (runtype == "bomd-electronic") {
      bomd.reset(new ElectronicBomd(qm_grd, mm_grd));
    }
    else if (runtype == "ehrenfest") {
      bomd.reset(new Ehrenfest(qm_grd, mm_grd));
    }
    else if (runtype == "diabatic_seam") {
      bomd.reset(new DiabaticSeam(qm_grd, mm_grd));
    }
    else {
      throw std::runtime_error("unknown runtype");
    }
    std::cerr << "Begining INAQS run of type " << runtype << std::endl;
    qmi = bomd->setup(fh, atomicnumbers, qm_crd, mm_crd, mm_chg);
    
    return bomd;
};

void
GifsImpl::setup_reader(ConfigBlockReader& reader) {
  //using types = ConfigBlockReader::types;
    //
    reader.add_entry("runtype", std::string("bomd"));
//    reader.add_entry("latoms", std::vector<int>{});
}


GifsImpl::GifsImpl(const char* file, size_t in_nqm, const int * ian, const double mass, const double length, const double time)
    : nqm{in_nqm}, nmm{0}, qm_atomicnumbers(in_nqm), qm_index(in_nqm), conv(mass, length, time)
{

  FileHandle fh{file};
  ConfigBlockReader reader;
  {
    std::vector<std::string> block_names = {
      "inaqs", // primary block first
      "gifs"   // legacy blocks follow
    };
    for (const auto & b: block_names){
      reader = ConfigBlockReader(b);
      setup_reader(reader);
      if (0 == reader.parse(fh)){
        if (b != block_names[0]){
          std::cerr << "DEPRECATION WARNING: you are not using the standard config block, ["
                    << block_names[0] << "]. Please update your inputs." << std::endl;
        }
        break;
      }
    }
  }

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
    //bomd = new BOMD(in_nqm, in_qmid);
    bomd = select_bomd(reader, fh,
                    qm_atomicnumbers,
                    qm_crd, 
                    mm_crd, 
                    mm_chg, 
                    qm_frc, 
                    mm_frc);
    // no link atoms
    arma::umat global_idx(0, 0);
    arma::vec factors(1);
    factors.resize(0);
    las.reset(LinkAtoms::with_const_factors(global_idx, factors));
};
//
void GifsImpl::destroy_instance() {
    if (instance_exists()) {
        delete impl;
    }
}

void Gifs::update_global_index(int* indexQM, int* indexMM) {
    impl->update_global_index(indexQM, indexMM);
}

void 
GifsImpl::update_global_index(int* indexQM, int* indexMM)
{
  //std::cout << "Update GLOBAL INDICES\n";
    for (arma::uword i=0; i<nqm; ++i) {
        qm_index[i] = indexQM[i];
        //std::cout << indexQM[i] << " ";
    };
    //std::cout << "\n";
    for (arma::uword i=0; i<nmm; ++i) {
        mm_index[i] = indexMM[i];
        //std::cout << indexMM[i] << " ";
    };
    //std::cout << "\n";
}
//
GifsImpl* GifsImpl::impl = nullptr;
//
template void GifsImpl::rescale_velocities(double total_energy, double* total_gradient, double* masses, double* velocities);
template void GifsImpl::rescale_velocities(float total_energy, float* total_gradient, float* masses, float* velocities);
template float GifsImpl::update_gradient(const float* in_qm_crd, const size_t* local_index, size_t in_nmm, const float* in_mm_crd, const float* in_mm_chg, float* in_qm_frc, float* in_mm_frc);
template double GifsImpl::update_gradient(const double* in_qm_crd, const size_t* local_index, size_t in_nmm, const double* in_mm_crd, const double* in_mm_chg, double* in_qm_frc, double* in_mm_frc);

  
