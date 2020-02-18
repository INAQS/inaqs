#include "fssh.hpp"
#include <armadillo>
#include <complex>

void FSSH::main(void){
  PropMap props{};
  props.emplace(QMProperty::qmgradient, {active_state}, &qm_grd);
  props.emplace(QMProperty::mmgradient, {active_state}, &mm_grd);
  props.emplace(QMProperty::energies,   {active_state}, &energy);
  props.emplace(QMProperty::wfoverlap,  &U);
  qm->get_properties(props);

  // FIXME: need to orthogonalize U first...
  
  //T = arma::logmat(U) / dtc;

  // compute min. time step for electronic propagation (eqs. 20, 21)
  double dtq_ = std::min(dtc,
			 std::min(0.02 / T.max(),
				  0.02 / (V.max() - arma::mean(V))));
  dtq = dtc / std::round(dtc / dtq_);
  const size_t n_steps = (size_t) dtc / dtq;

  
  // propagate electronic coefficients (for each set of amplitudes)
  for (size_t nt = 0; nt < n_steps; nt++){
    // propagate

    // see if we need to hop
    if (not hopping){
      arma::uword a = active_state;
      arma::vec g = -2 * arma::real(arma::conj(c) * c[a] * T.col(a)) * dtq / std::norm(c[a]);
      // set negative elements to 0
      //g.elem( find(g < 0) ).zeros();      
    }

  
  }
  
}
