#include "fssh.hpp"
#include "electronic.hpp"
#include <armadillo>
#include <complex>

double FSSH::gen_rand(void){
  static std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
  return uniform_distribution(mt64_generator);
}

// FIXME: need to parse our config
FSSH::FSSH(int nqm, const int * qmid, size_t min_state, size_t excited_states, size_t active_state, double dtc):
  BOMD(nqm, qmid), dtc {dtc}, min_state {min_state}, excited_states {excited_states}, active_state {active_state}
{
  int nstates = excited_states + 1 - min_state;
  if (nstates < 2){
    throw std::logic_error("Cannot run FSSH on a single surface!");
  }

  if ((active_state < min_state) || (excited_states < active_state)){
    throw std::range_error("Active state not in the range of computed states!");
  }
  
  // FIXME: some objects have excited + ground states, others have
  // excited + ground - minimum; this will screw up indexing.
  U.set_size(nstates);
  energy.set_size(excited_states + 1);
  T.set_size(nstates);
  V.set_size(nstates);

  // FIXME: want to configure multiple states
  // FIXME: should be able to configure (multiple) state(s)
  arma::cx_vec c_ (nstates, arma::fill::zeros);
  c_(active_state) = 1.0;
  c.set(c_);
  // c.set_size(nstates, 1);
  // c.zeros();
  // c(active_state,0) = 1.0;
  
  std::random_device rd;
  mt64_generator.seed(rd()); //FIXME: should be able to pass random seed
}; 

// For use in the update_gradient() call; Jain step 4
void FSSH::electonic_evolution(void){
  if (hopping){
    throw std::logic_error("Should never be hopping at start of electronic evolution!");
  }
  
  PropMap props{};
  props.emplace(QMProperty::qmgradient, {active_state}, &qm_grd);
  props.emplace(QMProperty::mmgradient, {active_state}, &mm_grd);
  props.emplace(QMProperty::energies,   {excited_states}, &energy); // R.B.: with any idx, energies will get *all* states
  props.emplace(QMProperty::wfoverlap,  &U);
  qm->get_properties(props);

  // FIXME: need to match previous phases and orthogonalize U first
  // We could monitor the magnitude of the imaginary component
  T = real(arma::logmat(U)) / dtc;
    
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
    // propagate all states simutaneously
    //c = electronic::rk4_advance(V - I*T, c, dtq);
    c.advance(V - I*T, dtq);
    
    // check for a hop unless we've already had one
    if (! hopping){
      const arma::uword a = active_state;

      // transition probabilities; eq. 12 has the conjugate on the wrong element; see Tully (1990). 
      arma::vec g = -2 * arma::real(arma::conj(c.col(0)) * c(a,0) * T.col(a)) * dtq / std::norm(c(a,0));
      // set negative elements to 0
      g.elem( find(g < 0) ).zeros();
      
      // index of largest transition probability
      arma::uword j = g.index_max(); 
      if (g[j] > gen_rand()){
	// will update these in the velocity_rescale call
	hopping = true;
	target_state = j;
      }
    }
  }
}

// For use within the velocity_rescale() call; Jain steps 5 & 6
void FSSH::attempt_hop(void){
  // just in-case this gets called when it shouldn't
  if (! hopping){
    return;
  }

  
  
  
  hopping = false;
}
