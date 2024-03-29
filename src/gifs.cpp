#include <vector>
#include "gifs_implementation.hpp"
//#include "gifs.hpp"

//FIXME: do this in a less terrible way xD
bool is_init = false;

#ifdef __cplusplus
extern "C" {
#endif 
//Don't inline these; they need to be exported and available to Gromacs etc...
void create_qm_interface(const char* file, int mdStepStart, double classicalTimeStep, size_t nqm, const int* qm_atomids)
{
  const double mass_unit = 1.660538921e-27; 
  const double length_unit = 1e-9; 
  const double time_unit = 1e-12;
  Gifs gifs_handle(file, mdStepStart, classicalTimeStep, nqm, qm_atomids, mass_unit, length_unit, time_unit);
  is_init = true;
}

  bool gifs_interface_is_ready(void){ return is_init; }

  float gifs_get_forces_float(const float* qm_crd,
                              const size_t* link_ids,
                              size_t nmm,
                              const float* mm_crd,
                              const float* mm_chg,
                              float* f_qm, float* f_mm)
  {
    Gifs gifs_handle;
    return gifs_handle.update_gradient(qm_crd,
                                       link_ids, nmm,
                                       mm_crd, mm_chg,
                                       f_qm, f_mm);
  };

  double gifs_get_forces_double(const double* qm_crd,
                                const size_t* link_ids,
                                size_t nmm,
                                const double* mm_crd,
                                const double* mm_chg,
                                double* f_qm, double* f_mm)
  {
    Gifs gifs_handle;
    return gifs_handle.update_gradient(qm_crd,
                                       link_ids, nmm,
                                       mm_crd, mm_chg,
                                       f_qm, f_mm);
  };
  
void
gifs_update_global_index(int* indexQM, int* indexMM) {
  Gifs gifs_handle;
  gifs_handle.update_global_index(indexQM, indexMM);
}

  void
  gifs_rescale_velocities_float(float total_energy, float* total_gradient, float* masses, float* velocities) {
    if (is_init){
      Gifs gifs_handle;
      gifs_handle.rescale_velocities(total_energy, total_gradient, masses, velocities);
    }
  }

  void
  gifs_rescale_velocities_double(double total_energy, double* total_gradient, double* masses, double* velocities) {
    if (is_init){
      Gifs gifs_handle;
      gifs_handle.rescale_velocities(total_energy, total_gradient, masses, velocities);
    }
  }

  
#ifdef __cplusplus
}

std::shared_ptr<QMInterface> gifs_QMInterface(void){
  Gifs gifs_handle;
  return gifs_handle.get_QMInterface();
}

void gifs_update_coords(const arma::mat & R){
  Gifs gifs_handle;
  gifs_handle.update_coords(R);
}
#endif 
