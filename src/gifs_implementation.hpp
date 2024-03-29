#ifndef GIFS_SH_GIFS_IMPLEMENTATION_HPP
#define GIFS_SH_GIFS_IMPLEMENTATION_HPP

#include <vector>
#include <memory>
#include <armadillo>
//
#include "bomd.hpp"
#include "conversion.hpp"
#include "linkatoms.hpp"
#include "configreader.hpp"
#include "qm_interface.hpp"
#include "inaqs_shared.hpp"

/* Virtual base class, defining all operations that 
 * should be called from the MD software */

// factory to create new bomd instance, based on an input file!


// consider using a reference to a method-static object in a static method like:
// https://codereview.stackexchange.com/questions/13059/thread-safe-singleton-class-using-stdshared-ptr-in-c11/13082#13082

class GifsImpl
{
public: 
    // creation
    static GifsImpl* get_instance(const char* file, int mdStepStart, double classicalTimeStep, size_t nqm, const int * qmid,
            const double mass, const double length, const double time);
    static GifsImpl* get_instance();
    // methods
    template<typename T>
    T update_gradient(const T* in_qm_crd, const size_t* local_index, 
		      size_t in_nmm, const T* in_mm_crd, const T* in_mm_chg, 
		      T* in_qm_frc, T* in_mm_frc);
    //
    template<typename T> 
    void rescale_velocities(T total_energy, T* total_gradient, T* masses, T* velocities);
    //
    void update_global_index(int* indexQM, int* indexMM);
    //
    static inline bool instance_exists() noexcept {
        if (impl == nullptr) {
            return false;
        }
        return true;
    }

  std::shared_ptr<QMInterface> get_QMInterface(void){return qmi;}
  void update_coords(const arma::mat & R);
  
private:
  GifsImpl(const char* file, int mdStepStart, double classicalTimeStep, size_t nqm, const int * qmids,
             const double mass, const double length, const double time);
    // private destructor
    static void destroy_instance(void);
    // Singleton Pointer
    static GifsImpl* impl; 
    //
    void setup_reader(ConfigBlockReader&);
    std::unique_ptr<BOMD> select_bomd(ConfigBlockReader& reader, FileHandle& fh,
                                      arma::uvec& atomicnumbers,
                                      arma::mat& qm_crd, arma::mat& mm_crd,
                                      arma::vec& mm_chg,
                                      arma::mat& qm_grd, arma::mat& mm_grd);

    // Local data!
    arma::uword nqm;
    arma::uword nmm;
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
    // velocity rescaling! 
    // all these are only "real" atoms: 
    // nqm + nmm without linkatom!
    arma::mat veloc{};
    arma::vec masses{};
    arma::mat total_gradient{};
    // Classes
    Conversion conv;

    std::unique_ptr<BOMD> bomd;
    std::unique_ptr<LinkAtoms> las;
    std::shared_ptr<QMInterface> qmi;
    std::shared_ptr<INAQSShared> shared;
};


class Gifs
{
public:
    // create instance, can only be called once!
    explicit Gifs(const char* file, int mdStepStart, double classicalTimeStep, size_t nqm, const int * qmid, const double mass, const double length, const double time) {
        if (!GifsImpl::instance_exists()) {
            impl = GifsImpl::get_instance(file, mdStepStart, classicalTimeStep, nqm, qmid, mass, length, time);
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
        return impl->update_gradient(qm_crd, local_index, 
                              nmm, mm_crd, mm_chg, 
                              qm_frc, mm_frc);
    };
    //
    template<typename T> 
    inline
    void rescale_velocities(T total_energy, T* total_gradient, T* masses, T* velocities) {
      impl->rescale_velocities(total_energy, total_gradient, masses, velocities);
    }
    //
    void update_global_index(int* indexQM, int* indexMM);

    std::shared_ptr<QMInterface> get_QMInterface(void){return impl->get_QMInterface();}
    void update_coords(const arma::mat & R){impl->update_coords(R);}

private:
    GifsImpl* impl{nullptr}; 
};

#endif // GIFS_SH_GIFS_CORE_H
