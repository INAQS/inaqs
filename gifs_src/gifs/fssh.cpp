#include "fssh.hpp"
#include "electronic.hpp"
#include "constants.hpp"
#include <armadillo>
#include <complex>
#include <cmath>
#include <iostream>

ConfigBlockReader
FSSH::setup_reader()
{
    using types = ConfigBlockReader::types;
    ConfigBlockReader reader{"fssh"};
    reader.add_entry("dtc", types::DOUBLE);
    reader.add_entry("delta_e_tol", 1e-4);

    reader.add_entry("amplitude_file", "cs.dat");

    // FIXME: make ConfigBlockReader complain if adding the same key twice!
    
    reader.add_entry("active_state", types::ULINT);

    std::random_device rd; // generate default random seed
    reader.add_entry("random_seed", (size_t) rd());
    return reader;
}


void
FSSH::get_reader_data(ConfigBlockReader& reader) {
  {
    double in_dtc;
    reader.get_data("dtc", in_dtc);
    // FIXME?: Should input time be in fs as it is (because our
    // interface is general) or in ns to match GMX?
    dtc = in_dtc * (1e-15 / AU2SI_TIME); // fs -> a.u.
  }

  reader.get_data("delta_e_tol", delta_e_tol);
  reader.get_data("active_state", active_state);
  
  /* added in BOMD::add_qm_keys() */
  reader.get_data("min_state", min_state);
  reader.get_data("excited_states", excited_states);

  reader.get_data("amplitude_file", amplitude_file);
  
  {
    size_t seed = 0; // a default seed is set in setup_reader(); this value will not be propagated.

    // FIXME: ConfigReader won't throw an error if the wrong type is passed in
    reader.get_data("random_seed", seed);
    mt64_generator.seed(seed);
    arma::arma_rng::set_seed(seed);

    std::cerr << "FSSH random seed=" << seed << std::endl;
  }

  shstates = excited_states + 1 - min_state;
  if (shstates < 2){
    throw std::logic_error("Cannot run FSSH on a single surface!");
  }

  if ((active_state < min_state) || (excited_states < active_state)){
    throw std::range_error("Active state not in the range of computed states!");
  }

  active_state -= min_state;

  if (min_state != 1){
    std::cerr << "min_state=" << min_state << std::endl;
    throw std::runtime_error("minimum states other than 1 not supported!");
  }

  /*
    R.B. energy has excited_states + 1 (ground) states, while all
    other objects have excited_states + 1 - min_state (shstates); this
    can screw up indexing.

      active_state \on [0, shstates]

    can index U, T, V, c but cannot index energy which has shstates +
    min_state elements.  Indexing properly by taking:

      c(i) -> energy(i+min_state).

    Recall also that all requests for state properties to Q-Chem (or
    any other QM_Interface) require i+min_state as well.
  */

  energy.set_size(excited_states + 1);

  U.set_size(shstates, shstates);
  T.set_size(shstates, shstates);
  V.set_size(shstates, shstates);

  /*
    FIXME: should be able to configure (multiple) initial electronic
    state(s). Can use DVEC?
  */
  c.reset(shstates, 1, active_state);
};


FSSH::FSSH(arma::mat& qm_grd, arma::mat& mm_grd):
  BOMD{qm_grd, mm_grd}
{};


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
  if (p.has_nan()){throw std::logic_error("P is malformed; it contains NaN!");}
  if (arma::any(p < 0)){throw std::logic_error("P cannot have negative elements!");}
  if (std::abs(arma::sum(p) - 1) > 1e-8){ throw std::logic_error("P is not normed!");}

  const double zeta = arma::randu();
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

  arma::mat V = diagmat(energy.tail(shstates));

  /*
    Since the off-diagonal elements of V are always 0 in this
    implementation, we work only with the diagonal of V. Equation 20
    seems to have been developed for a system with positive energy
    eigenvalues. A constant shift of energy scale to the negative (Vii
    = Vii - max(V.diag) can cause Eq. 20 to fail.

    FIXME: construct an alternate to Eq. 20 that is applicable to Vii < 0
  */

  // compute max time step for electronic propagation (eqs. 20, 21)
  {
    double dtq_ = std::min(dtc,
			   std::min(0.02 / T.max(),
				    0.02 / arma::max( V.diag() - arma::mean(V.diag()) )
				    )
			   );

    dtq = dtc / std::round(dtc / dtq_);
  }
  const size_t n_steps = (size_t) dtc / dtq;
  const std::complex<double> I(0,1);

  // Propagate electronic coefficients and compute hopping probabilities for dtc
  for (size_t nt = 0; nt < n_steps; nt++){
    /*
      Propagate all states simultaneously. Recall that
      Electronic::advance(H, dt) uses rk4 to propagate the internally
      held coefficients, c, according the Hamiltonian H for time dt
    */

    // FIXME: should use rk4
    //c.advance_rk4(V - I*T, dtq);
    c.advance_exact(V - I*T, dtq);

    // Check for a hop unless we've already had one
    if (! hopping){
      const arma::uword a = active_state;

      /*
	Compute transition probabilities. Jain (2016) eq. 12 has the
        conjugate on the wrong element; cf. Tully (1990) -- discussion
        with Zeyu Zhou
      */
      arma::vec g = -2 * arma::real(c(a) * arma::conj(c()) % T.col(a)) * dtq / std::norm(c(a));
      // set negative elements to 0
      g.elem( arma::find(g < 0) ).zeros();

      /*
	Ensure that g is normed by adding any residual density to the
	active state. This maintains the correct transition
	probability to all states.
      */
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
void FSSH::hop_and_scale(arma::mat &velocities, arma::vec &mass){
  if (! hopping){
    throw std::logic_error("Should not be attempting a hop right now!");
  }

  /*
    Recall that in our implementation, the kinetic energy reservoir to
    balance a hop includes the MM atoms, a region of arbitrary
    size. However, the kinetic energy is only available in proportion
    to the NAC vector, which is properly computed over the MM region
    (Thank you, Ou Qi!)  and, we believe, decays with distance. See
    the calculation of vd below.
  */

  if (mass.n_elem != NQM() + NMM()){
    throw std::range_error("mass of improper size!");
  }

  nac.set_size(3, NQM() + NMM());

  PropMap props{};
  props.emplace(QMProperty::nacvector, {min_state + active_state, min_state + target_state}, &nac);

  qm->get_properties(props);

  // Make 3N vector versions of the NAC, velocity, and mass
  // FIXME: How do each of these interact with link atoms?
  const arma::vec nacv(nac.memptr(), 3 * (NQM() + NMM()), false, true);
  arma::vec vel(velocities.memptr(), 3 * (NQM() + NMM()), false, true);


  arma::vec m (3 * (NQM() + NMM()), arma::fill::zeros);
  for(arma::uword i = 0 ; i < m.n_elem ; i++){
    m[i] =  mass[i/3]; // without bounds check
    //m(i) =  mass(i/3); // with bounds check
  }

  // Unlike all other state-properties, must use min_state as floor for indexing into energy
  double deltaE = energy(min_state + target_state) - energy(min_state + active_state);

  // The hopping energy-conservation equations are documented in documented Vale's GQSH notes Sept 11, 2020

  double dmv = arma::as_scalar(nacv.t() * (vel % m));
  double dmd = arma::as_scalar(nacv.t() * (nacv % m));

  double discriminant = (dmv/dmd)*(dmv/dmd) - 2*deltaE/dmd;
  if (discriminant > 0){
    std::cerr << "Hop, " << active_state + min_state << "->" << target_state + min_state << " succeeds; energy difference = " << deltaE << "."<< std::endl;
    // test the sign of dmv to pick the root yielding the smallest value of alpha
    double alpha = (dmv > 0 ? 1.0 : -1.0) * std::sqrt(discriminant) - (dmv/dmd);
    vel = vel + alpha * nacv;
    active_state = target_state;
  }
  else{  // frustrated hop
    std::cerr << "Hop is frustrated---will remain on " << active_state + min_state << "; ";
    // compute gradient of target_state to see if we reverse
    arma::mat qmg_new, mmg_new;
    qmg_new.set_size(3,NQM());

    /*
      There is, perhaps, a modest performance advantage to be gained
      by saving this gradient for use below in
      backpropagate_gradient_velocities(). This would be most important for
      situations where we had many frustrated hops. The infrastructure
      to track which gradients are up-to-date is a complication for
      later.
    */

    props = {};
    props.emplace(QMProperty::qmgradient, {min_state + target_state}, &qmg_new);
    if (NMM() > 0){
      mmg_new.set_size(3, NMM());
      props.emplace(QMProperty::mmgradient, {min_state + target_state}, &mmg_new);
    }
    qm->get_properties(props);

    // 3N vector version of new gradient
    arma::vec gradv(3 * (NQM() + NMM()));
    gradv.head(3*NQM()) = arma::vectorise(qmg_new);
    if (NMM() > 0){
      gradv.tail(3*NMM()) = arma::vectorise(mmg_new);
    }

    /*
      Velocity reversal along nac as per Jasper, A. W.; Truhlar,
      D. G. Chem. Phys. Lett. 2003, 369, 60--67 c.f. eqns. 1 & 2

      In Jain (2016) a second criterion was imposed. But, in January
      2020 A. Jain indicated to JES that this was not necessary. We
      follow the original Jasper-Truhlar prescription in line with
      Jain's updated advice.
    */
    if (arma::as_scalar((-gradv.t() * nacv)*(vel.t() * nacv)) < 0){
      std::cerr << "velocities reversed." << std::endl;
      arma::vec nacu = arma::normalise(nacv);
      vel = vel - 2.0 * nacu * nacu.t() * vel;
    }
    else{
      std::cerr << "velocities, unchanged." << std::endl;
      // Ignore the unsuccessful hop; active_state remains unchanged
    }
  }
  hopping = false;
}



/*
  Call this function when energy fluctuation tolerance in the MD
  driver is exceeded and we've had a hop. Update the current surface
  and to the total gradient, add (new-old). Similarly, take a step in
  velocity space backwards along the old gradient and then forwards
  along the new one.

  N.B.: When working with Gromacs, it is not necessary to do any
  back-propagation in positions, only velocities. Gromacs's
  velocity-Verlet integrator splits velocity integration over 2 steps
  (so that all of the tooling for leap-frog still works) and then does
  position integration in a single step afterwards. See
  gifs/docs/GromacsVVImplementation.jpg for notes.
*/
void FSSH::backpropagate_gradient_velocities(arma::mat &total_gradient, arma::mat &velocities, arma::vec &masses){
  PropMap props = {};
  arma::mat qmg_new, mmg_new;
  qmg_new.set_size(3,NQM());

  props.emplace(QMProperty::qmgradient, {min_state + target_state}, &qmg_new);
  if (NMM() > 0){
    mmg_new.set_size(3, NMM());
    props.emplace(QMProperty::mmgradient, {min_state + target_state}, &mmg_new);
  }
  qm->get_properties(props);

  // 3xN  mass matrix
  arma::mat m (3, (NQM() + NMM()), arma::fill::zeros);
  for (arma::uword i = 0; i < masses.n_elem; i++){
    m.col(i) = masses(i);
  }

  auto rows=arma::span(arma::span::all); auto cols = arma::span(0, NQM() - 1);
  total_gradient(rows, cols) += qmg_new - qm_grd;
  // a = -g/m
  velocities(rows, cols) += -1.0*(qmg_new - qm_grd)/m(rows, cols)*dtc/2;
  qm_grd = qmg_new;

  if (NMM() > 0){
    cols = arma::span(NQM(), NQM() + NMM() - 1);
    total_gradient(rows, cols) += mmg_new - mm_grd;
    velocities(rows, cols) += -1.0*(mmg_new - mm_grd)/m(rows, cols)*dtc/2;
    mm_grd = mmg_new;
  }
}


// This is our primary hook into the Gromacs (or other) MD loop
double FSSH::update_gradient(void){
  qm->update();

  PropMap props{};
  props.emplace(QMProperty::qmgradient, {min_state + active_state}, &qm_grd);
  props.emplace(QMProperty::mmgradient, {min_state + active_state}, &mm_grd);
  // R.B.: with any idx, energies will get all (excited_states + 1) states
  props.emplace(QMProperty::energies,   {excited_states}, &energy);

  qm->get_properties(props);

  std::ofstream output(amplitude_file, std::ios_base::app);
  c().t().print(output);
  output.close();

  electonic_evolution();

  return energy(min_state + active_state);
}


// This hook comes after the first MD half-step (and before constraint forces are calculated in gromacs)
bool FSSH::rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double total_energy){
  BOMD::rescale_velocities(velocities, masses, total_gradient, total_energy); // call parent to update edrift
  if (!hopping){
    // nothing to do
    return false;
  }
  else{
    /*
      FIXME: For QM/MM some species hidden from the QM, this energy is
      for the TOTAL system, which is probably not what we want.
    */
    if(std::abs(edrift) > delta_e_tol){
      std::cerr <<  "Trivial crossing: " << active_state + min_state << "->" << target_state + min_state << std::endl;
      // Need to update global gradient & velocities
      backpropagate_gradient_velocities(total_gradient, velocities, masses);
      active_state = target_state;
      hopping = false;
    }
    else{
      std::cerr <<  "Attempting hop: " << active_state + min_state << "->" << target_state + min_state << std::endl;
      hop_and_scale(velocities, masses);
      hopping = false;
    }
  }
  return true;
}
