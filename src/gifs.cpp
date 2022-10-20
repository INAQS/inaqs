#include <vector>
#include "gifs_implementation.hpp"
//#include "gifs.hpp"

//FIXME: do this in a less terrible way xD
bool is_init = false;

#ifdef __cplusplus
extern "C" {
#endif 
//Don't inline these; they need to be exported and available to Gromacs etc...
void create_qm_interface(const char* file, size_t nqm, const int* qm_atomids)
{
  const double mass_unit = 1.660538921e-27; 
  const double length_unit = 1e-9; 
  const double time_unit = 1e-12;
  Gifs gifs_handle(file, nqm, qm_atomids, mass_unit, length_unit, time_unit);
  is_init = true;
}


float gifs_get_forces(const float* qm_crd, 
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

void
gifs_update_global_index(int* indexQM, int* indexMM) {
  Gifs gifs_handle;
  gifs_handle.update_global_index(indexQM, indexMM);
}

void
gifs_rescale_velocities(float total_energy, float* total_gradient, float* masses, float* velocities) {
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
#endif 
