#include "fssh.hpp"
#include "electronic.hpp"
#include <armadillo>
#include <complex>
#include <cmath>

double FSSH::gen_rand(void){
  static std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
  return uniform_distribution(mt64_generator);
}


// FIXME: need to parse our config
FSSH::FSSH(int nqm, const int * qmid, size_t min_state, size_t excited_states, size_t active_state, double dtc):
  BOMD(nqm, qmid), dtc {dtc},
  min_state {min_state}, excited_states {excited_states}, active_state {active_state},
  c {excited_states + 1 - min_state}
{
  int nstates = excited_states + 1 - min_state;
  if (nstates < 2){
    throw std::logic_error("Cannot run FSSH on a single surface!");
  }

  if ((active_state < min_state) || (excited_states < active_state)){
    throw std::range_error("Active state not in the range of computed states!");
  }

  // FIXME: need to pass min_state to QM_Interface

  // FIXME: some objects have excited + ground states, others have
  // excited + ground - minimum; this will screw up indexing.
  U.set_size(nstates);
  energy.set_size(excited_states + 1);
  T.set_size(nstates);
  V.set_size(nstates);

  // FIXME: should be able to configure (multiple) state(s)
  // R.B.: c is of type Electronic, not an arma object
  arma::cx_vec c_ (nstates, arma::fill::zeros);
  c_(active_state) = 1.0;
  c.set(c_);
  
  std::random_device rd;
  mt64_generator.seed(rd()); //FIXME: should be able to pass random seed
  //FIXME: attempting to set the armadillo random seed errors: ‘arma_rng’ has not been declared
  //arma_rng::set_seed_random(); //arma_rng::set_seed(rd());
}; 


/*
  For a discrete probability distribution Pi = p(i), returns the index
  of an element randomly selected by sampling the distribution.

  Examples:
  * sample_discrete({0.5, 0.5}) returns 0 or 1 with probability 0.5
  * sample_discrete({0.0, 0.0, 1.0}) returns 2 with probability 1.0

  Note:
  While there does exist std::discrete_distribution, I think that the
  below implementation is superior in that it does not require
  instantiating a new distribution for every round.
*/
arma::uword FSSH::sample_discrete(const arma::vec &p){
  /*
    Sanity checks; require:
    * Pi >= 0 forall i 
    * Sum(Pi) == 1
  */
  if (arma::any(p < 0)){throw std::logic_error("P cannot have negative elements!");}
  if (std::abs(arma::sum(p) - 1) > 1e-8){ throw std::logic_error("P is not normed!");}

  const double zeta = gen_rand();
  return as_scalar(arma::find(arma::cumsum(p) > zeta, 1));
}


// For use in the update_gradient() call; Jain step 4
void FSSH::electonic_evolution(void){
  if (hopping){
    throw std::logic_error("Should never be hopping at start of electronic evolution!");
  }
  
  PropMap props{};
  props.emplace(QMProperty::wfoverlap,  &U);
  qm->get_properties(props);

  Electronic::phase_match(U);
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
  for (size_t nt = 0; nt < n_steps; nt++){
    // propagate all states simutaneously
    c.advance(V - I*T, dtq);

    // check for a hop unless we've already had one
    if (! hopping){
      const arma::uword a = active_state;

      // transition probabilities; eq. 12 has the conjugate on the wrong element; see Tully (1990).
      arma::vec g = -2 * arma::real(arma::conj(c()) * c(a) * T.col(a)) * dtq / std::norm(c(a));
      // set negative elements to 0
      g.elem( arma::find(g < 0) ).zeros();

      // FIXME: is this norming the correct thing to do? Should we rather be monitoring the norm of g?
      // ensure that g is normed by adding any residual desnity to the active state.
      g(a) += 1.0 - arma::sum(g);

      // randomly select an element from the discrete distribution represented by g
      arma::uword j = sample_discrete(g);
      if (a != j){
      	// will update these in the velocity_rescale call
      	hopping = true;
      	target_state = j;
      }
    }
  }
}


// For use within the velocity_rescale() call; Jain steps 5 & 6
void FSSH::hop_and_scale(arma::vec &vel, arma::vec &mass){
  if (! hopping){
    throw std::logic_error("Should not be attempting a hop right now!");
  }

  // FIXME: MFSJM: does this accord with your plan to deal with link atoms?
  // Perhaps these should be members of BOMD?
  arma::uword NQM = qm_grd.n_cols;
  arma::uword NMM = mm_grd.n_cols;

  if (mass.n_elem != NQM + NMM){
    throw std::range_error("mass of improper size!");
  }

  nac.set_size(3, NQM + NMM);
  
  PropMap props{};
  //FIXME: Make sure you don't need to add min_state to the below
  props.emplace(QMProperty::nacvector, {active_state, target_state}, &nac);
  
  qm->get_properties(props);

  // Make 3N vector versions of the NAC and mass
  const arma::vec nacv(nac.memptr(), 3 * (NQM + NMM), false, true);

  arma::vec m (3 * (NQM + NMM), arma::fill::zeros);
  for (arma::uword i = 0; i < mass.n_elem; i++){
    m.subvec(3 * i, 3*(i+1)) = mass(i);
  }


  // Unlike all other state-properties, must use min_state as floor for indexing into energy
  double deltaE = energy(min_state + target_state) - energy(min_state + active_state);
  
  double vmd = arma::as_scalar((vel % m)  * nacv.t());
  double dmd = arma::as_scalar((nacv % m) * nacv.t());

  double discriminant = (vmd/dmd)*(vmd/dmd) - 2*deltaE/dmd;
  if (discriminant > 0){  // hop succeeds
    // test the sign of vmd to pick the root yeilding the smallest value of alpha
    double alpha = (vmd > 0 ? 1.0 : -1.0) * std::sqrt(discriminant) - (vmd/dmd);
    // FIXME: verify the dimension (units) of the NAC as calculated by qchem; do we compute NAC or DC?
    vel = vel + alpha * nacv;
    active_state = target_state;
  }
  else{  // frustrated hop
    // compute gradient of target_state to see if we reverse
    arma::mat qmg_new, mmg_new;
    qmg_new.set_size(3,NQM);
  
    props = {};
    props.emplace(QMProperty::qmgradient, {target_state}, &qmg_new);
    if (NMM > 0){
      mmg_new.set_size(3, NMM);
      props.emplace(QMProperty::mmgradient, {target_state}, &mmg_new);
    }
    qm->get_properties(props);

    // 3N vector version of gradient
    arma::vec gradv(3 * (NQM + NMM));
    gradv.subvec(0, 3*NQM) = arma::vectorise(qmg_new);
    if (NMM > 0){
      gradv.subvec(3*NQM + 1, 3 * (NQM + NMM)) = arma::vectorise(mmg_new);
    }
    
    /*
      Velocity reversal along nac as per Jasper, A. W.; Truhlar,
      D. G. Chem. Phys. Lett. 2003, 369, 60--67 c.f. eqns. 1 & 2

      In the original Jain paper a second criterion was imposed. But,
      in January 2020 A. Jain indicated to Joe that this was not
      necesary. We follow the original Jasper-Truhlar perscription in
      line with Jain's updated advice.
    */
    if (arma::as_scalar((-gradv * nacv)*(vel * nacv)) < 0){
      arma::vec nacu = arma::normalise(nacv);
      vel = vel - 2.0 * nacu * nacu.t() * vel;
    }
    else{
      // Ignore the unsuccessful hop; active_state remains unchanged 
    }
  }
  hopping = false;
}



/*
  Call this function when energy fluctuation tolerance in the MD
  driver is exceeded and we've had a hop. Update the current surface
  and to the total gradient, add (new-old).
*/
void FSSH::update_MD_global_gradient(void){
  PropMap props = {};
  arma::mat qmg_new, mmg_new;
  qmg_new.set_size(3,NQM);
  

  props.emplace(QMProperty::qmgradient, {target_state}, &qmg_new);
  if (NMM > 0){
    mmg_new.set_size(3, NMM);
    props.emplace(QMProperty::mmgradient, {target_state}, &mmg_new);
  }
  qm->get_properties(props);

  // add this to whatever representation of the total gradient that we have
  //qmg_new - qm_grd;
  //mmg_new - mm_grd;

  //FIXME: make sure this does a copy
  qm_grd = qmg_new;
  mm_grd = mmg_new;
}


double FSSH::update_gradient(void){
  //FIXME: should figure out how to call qm->update(void) automatically
  //qm->update();

  PropMap props{};
  props.emplace(QMProperty::qmgradient, {active_state}, &qm_grd);
  props.emplace(QMProperty::mmgradient, {active_state}, &mm_grd);
  props.emplace(QMProperty::energies,   {excited_states}, &energy); // R.B.: with any idx, energies will get *all* states

  qm->get_properties(props);
  electonic_evolution();

  return energy(min_state + active_state);
}


void FSSH::velocityrescale(void){  
  if (hopping){/*
    if(delta_e_from_md > md_tolerance){
      // trivial crossing; need to update global gradient
      update_md_global_gradient();
      active_state = target_state;
      hopping = false;
    }
    else{
      hop_and_scale(arma::vec &vel, arma::vec &mass);
      hopping = false;
    }
  */
  }
  else{
    // nothing to do
  }
}
