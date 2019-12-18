#include "qm.hpp"

class BOMD
{
public:
  BOMD(std::string config);

  /*
    for each of these templated functions, will need to:
    static_assert<std::is_floating_point<T>::value, "Error Msg">;
  */
  
  template<typename T>
  T get_gradient(T* qm_crd, T* mm_crd,
		 T* mm_chg, T* qm_gradient,
		 T* mm_gradient);

  template<typename T>
  T rescale_velocities(T* total_gradient, T* masses, T* velocities);

protected:
    QMInterface* QM;
    // fixed size
    std::vector<double> qm_grd;
    // flexible size
    std::vector<double> mm_grd;
};
