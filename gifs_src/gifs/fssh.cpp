#include "fssh.hpp"
#include "electronic.hpp"
#include <armadillo>
#include <complex>

double FSSH::gen_rand(void){
  static std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
  return uniform_distribution(mt64_generator);
}

FSSH::FSSH(int nqm, const int * qmid): // FIXME: need to parse our config
  BOMD(nqm, qmid){
  std::random_device rd;
  mt64_generator.seed(rd()); //FIXME: should be able to pass random seed
}; 

void FSSH::main(void){
  PropMap props{};
  props.emplace(QMProperty::qmgradient, {active_state}, &qm_grd);
  props.emplace(QMProperty::mmgradient, {active_state}, &mm_grd);
  props.emplace(QMProperty::energies,   {excited_states}, &energy);
  props.emplace(QMProperty::wfoverlap,  &U);
  qm->get_properties(props);

  // FIXME: need to orthogonalize U first?
  // FIXME: need to match previous phases
  arma::mat T_ = real(arma::logmat(U)) / dtc;
  arma::mat T = T_.submat(min_state, min_state, excited_states, excited_states);
    
  arma::mat V = diagmat(energy.subvec(min_state, excited_states));
  
  // compute min. time step for electronic propagation (eqs. 20, 21)
  double dtq_ = std::min(dtc,
			 std::min(0.02 / T.max(),
				  0.02 / (V.max() - arma::mean(V.diag()))));
  dtq = dtc / std::round(dtc / dtq_);
  const size_t n_steps = (size_t) dtc / dtq;

  const std::complex<double> I(0,1);
  
  // propagate electronic coefficients (for each set of amplitudes)
  hopping = false; // FIXME: reset in velocity_rescale
  for (size_t nt = 0; nt < n_steps; nt++){
    // propagate
    c = electronic::rk4_step(I*V - T, c, dtq);
    
    // check for a hop unless we've already had one
    if (! hopping){
      const arma::uword a = active_state;

      // transition probabilities; eq. 12 has the conjugate on the wrong element; see Tully (1990). 
      arma::vec g = -2 * arma::real(arma::conj(c) * c[a] * T.col(a)) * dtq / std::norm(c[a]);
      // set negative elements to 0
      g.elem( find(g < 0) ).zeros();
      
      // index of largest transition probability
      arma::uword j = g.index_max(); 
      if (g[j] > gen_rand()){
	hopping = true;
	target_state = j;
      }
    }

  }
  
}
