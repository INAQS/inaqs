#include "dynamics.hpp"

namespace BOMD{

  BOMD(std::string config){
  }

  template<typename T>
  T get_gradient(T* qm_crd, T* mm_crd,
		 T* mm_chg, T* qm_gradient,
		 T* mm_gradient){
    static_assert<std::is_floating_point<T>::value, "Error Msg">;
  }

  /*No velocity resacle for BOMD*/
  template<typename T> T rescale_velocities(T* total_gradient, T* masses, T* velocities){}

}

