#include "fssh.hpp"
#include "electronic.hpp"
#include "constants.hpp"
#include "decoherence_afssh.hpp"
#include "util.hpp"
#include <armadillo>
#include <complex>
#include <cmath>
#include <iostream>

ConfigBlockReader FSSH::setup_reader()
{
    using types = ConfigBlockReader::types;
    ConfigBlockReader reader{"fssh"};
    reader.add_entry("dtc", types::DOUBLE);
    reader.add_entry("delta_e_tol", 1e-4);

    reader.add_entry("amplitude_file", "cs.dat");
    reader.add_entry("decoherence", "");
    // FIXME: make ConfigBlockReader complain if adding the same key twice!

    std::random_device rd; // generate default random seed
    reader.add_entry("random_seed", (size_t) rd());
    {
      std::vector<std::complex<double>> cs {};
      reader.add_entry("amplitudes", cs);
    }
    return reader;
}


void FSSH::get_reader_data(ConfigBlockReader& reader) {
  {
    double in_dtc;
    reader.get_data("dtc", in_dtc);
    // FIXME?: Should input time be in fs as it is (because our
    // interface is general) or in ns to match GMX?
    dtc = in_dtc * (1e-15 / AU2SI_TIME); // fs -> a.u.
  }

  reader.get_data("delta_e_tol", delta_e_tol);
  reader.get_data("amplitude_file", amplitude_file);

  /* added in BOMD::add_qm_keys() */
  reader.get_data("min_state", min_state);


  /* 
     Must set rng seed before decoherence because Decoherence will use
     armadillo random number generator
  */
  {
    size_t seed = 0; // a default seed is set in setup_reader(); this value will not be propagated.

    // FIXME: ConfigReader won't throw an error if the wrong type is passed in
    reader.get_data("random_seed", seed);
    arma::arma_rng::set_seed(seed);

    std::cerr << "[FSSH] random_seed = " << seed << std::endl;
  }

  {
    size_t excited_states;
    /* added in BOMD::add_qm_keys() */
    reader.get_data("excited_states", excited_states);

    shstates = excited_states + 1 - min_state;
  }

  if (shstates < 2){
    throw std::logic_error("Cannot run FSSH on a single surface!");
  }

  if (!(min_state <= active_state && active_state <= min_state + shstates)){
    throw std::range_error("Active state not in the range of hopping states!");
  }

  {
    std::string decoherence_in;
    reader.get_data("decoherence", decoherence_in);

    if (decoherence_in == "") {
      // do nothing; already nullptr
    } else if (decoherence_in == "afssh" ||
               decoherence_in == "jain2016" ) {
      decoherence = new AFSSH(&qm, dtc, min_state, shstates, NQM(), NMM());
    }
    else{
      throw std::runtime_error("Decoherence class, '" + decoherence_in + "', not recognized!");
    }
  }

  active_state -= min_state;

  energy.set_size(shstates);

  U.set_size(shstates, shstates);
  T.set_size(shstates, shstates);
  V.set_size(shstates, shstates);

  /*
    FIXME: should be able to configure (multiple) initial electronic
    state(s). Can use DVEC?
  */

  {
    std::vector<std::complex<double>> cs_vec {};
    reader.get_data("amplitudes", cs_vec);
    if (cs_vec.size() > 0){
      if (cs_vec.size() != (arma::uword) shstates){
        throw std::runtime_error("Number of Amplitudes do not match number of hopping surfaces!");
      }
      double norm = 0;
      for (const auto& c: cs_vec){
        norm += std::norm(c);
      }
      std::cerr << "[FSSH] amplitude norm = " << norm << std::endl;
      if (std::abs(1.0 -norm) > 1e-6){
        throw std::runtime_error("Amplitudes are not normed; check your input!");
      }
      arma::cx_vec cs_arma(shstates);
      for (arma::uword i = 0; i < (arma::uword) shstates; i++){
        cs_arma(i) = cs_vec[i];
      }
      c = Electronic(cs_arma);
    }
    else{
      c.reset(shstates, 1, active_state);
    }
  }
}


// For use in the update_gradient() call; Jain step 4
void FSSH::electronic_evolution(void){
  if (hopping){
    throw std::logic_error("Should never be hopping at start of electronic evolution!");
  }

  PropMap props{};
  props.emplace(QMProperty::wfoverlap, &U);
  qm->get_properties(props);

  Electronic::phase_match(U);
  T = arma::real(arma::logmat(U)) / dtc;

  // R.B.: energy was updated in update_gradient()
  V = diagmat(energy);

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
      arma::uword j = util::sample_discrete(g);
      if (a != j){
      	// will update these in the velocity_rescale call
      	hopping = true;
      	target_state = j;
      }
    }  // end hopping check
  }  // end loop over dtq
}


// For use within the velocity_rescale() call; Jain steps 5 & 6
// returns energy gap if the hop succeeds without frustration
double FSSH::hop_and_scale(arma::mat &total_gradient, arma::mat &velocities, const arma::vec &m){
  if (! hopping){
    throw std::logic_error("Should not be attempting a hop right now!");
  }

  bool hop_succeeds = false;
  
  /*
    Recall that in our implementation, the kinetic energy reservoir to
    balance a hop includes the MM atoms, a region of arbitrary
    size. However, the kinetic energy is only available in proportion
    to the NAC vector, which is properly computed over the MM region
    (Thank you, Ou Qi!)  and, we believe, decays with distance. See
    the calculation of vd below.
  */

  nac.set_size(3, NQM() + NMM());

  arma::mat qmg_new, mmg_new;
  qmg_new.set_size(3,NQM());

  {
    PropMap props{};
    props.emplace(QMProperty::nacvector, {min_state + active_state, min_state + target_state}, &nac);
    props.emplace(QMProperty::qmgradient, {min_state + target_state}, &qmg_new);

    if (NMM() > 0){
      mmg_new.set_size(3, NMM());
      props.emplace(QMProperty::mmgradient, {min_state + target_state}, &mmg_new);
    }

    qm->get_properties(props);
  }
  
  // Make 3N vector versions of the NAC, velocity, and new
  // gradient. (m comes in as 3N.)
  // FIXME: How do each of these interact with link atoms?
  const arma::vec & nacv = nac.as_col();

  // FIXME: find a safer way to have a writeable velocity view
  arma::vec vel(velocities.memptr(), 3 * (NQM() + NMM()), false, true);

  arma::vec gradv(3 * (NQM() + NMM()));
  gradv.head(3*NQM()) = arma::vectorise(qmg_new);
  if (NMM() > 0){
    gradv.tail(3*NMM()) = arma::vectorise(mmg_new);
  }

  double deltaE = energy(target_state) - energy(active_state);

  /*
    The hopping energy-conservation equations are documented in
    Vale's GQSH notes dated May 14, 2021.
  */

  double vd  = arma::as_scalar(vel.t() * nacv);
  double dmd = arma::as_scalar(nacv.t() * (nacv / m));
  
  double discriminant = (vd/dmd)*(vd/dmd) - 2*deltaE/dmd;
  if (discriminant > 0){
    hop_succeeds = true;
    std::cerr << "[FSSH] Hop, " << active_state + min_state << "->"
              << target_state + min_state << " succeeds; energy difference = "
              << deltaE << std::endl;

    // test the sign of dmv to pick the root yielding the smallest value of alpha
    double alpha = (vd > 0 ? 1.0 : -1.0) * std::sqrt(discriminant) - (vd/dmd);
    vel = vel + alpha * (nacv / m); // recall the nacv has dimension of momentum
    active_state = target_state;
    
    // Update the gradient with the new surface so that GMX can take its second step
    auto rows=arma::span(arma::span::all); auto cols = arma::span(0, NQM() - 1);
    total_gradient(rows, cols) += qmg_new - qm_grd;
    qm_grd = qmg_new;

    if (NMM() > 0){
      cols = arma::span(NQM(), NQM() + NMM() - 1);
      total_gradient(rows, cols) += mmg_new - mm_grd;
      mm_grd = mmg_new;
    }
  }
  else{
    hop_succeeds = false;
    std::cerr << "[FSSH] Hop is frustrated---will remain on " << active_state + min_state << "; ";
    /*
      Momentum reversal along nac as per Jasper, A. W.; Truhlar,
      D. G. Chem. Phys. Lett. 2003, 369, 60--67 c.f. eqns. 1 & 2

      In Jain (2016) a second criterion was imposed. But, in January
      2020 A. Jain indicated to JES that this was not necessary. We
      follow the original Jasper-Truhlar prescription in line with
      Jain's updated advice.
    */
    if (arma::as_scalar((-gradv.t() * nacv)*(nacv.t() * (vel % m))) < 0){
      std::cerr << "velocities reversed." << std::endl;
      const arma::vec nacu = arma::normalise(nacv);
      //vel = vel - 2.0 * (nacu * nacu.t() * (vel % m))/m;
      vel = vel - 2.0 * (nacu / m) * (nacu.t() * (vel % m));
    }
    else{
      std::cerr << "velocities, unchanged." << std::endl;
      // Ignore the unsuccessful hop; active_state remains unchanged
    }
  }

  hopping = false;  // update class-level state
  return hop_succeeds ? deltaE : 0;
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

  Conversation with Amber in Marhc 2021 indicated that it was faster
  and didn't rely on an ad hoc selection of energy tolerance to simply
  converge properties in dtc directly rather than use this scheme.

  void FSSH::backpropagate_gradient_velocities(
    arma::mat &total_gradient,
    arma::mat &velocities,
    arma::vec &masses);
*/


// This is our primary hook into the Gromacs (or other) MD loop
double FSSH::update_gradient(void){
  qm->update();

  // get gradients and energies
  {
    PropMap props{};
    props.emplace(QMProperty::qmgradient, {min_state + active_state}, &qm_grd);
    props.emplace(QMProperty::mmgradient, {min_state + active_state}, &mm_grd);
    props.emplace(QMProperty::energies, util::range(min_state, min_state + shstates), &energy);
    qm->get_properties(props);
  }
  
  // write amplitudes
  {
    std::ofstream output(amplitude_file, std::ios_base::app);
    output << active_state + min_state << " ";
    c().st().print(output);
    output.close();
  }

  electronic_evolution();

  return energy(active_state);
}


/*
  rescale_velocities() comes after the first MD half-step (and before
  constraint forces are calculated in gromacs)
*/

// FIXME: should alter velocity rescale interface so we take inverse masses as 3N vector
bool FSSH::rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double total_energy){
  // call parent to update edrift
  BOMD::rescale_velocities(velocities, masses, total_gradient, total_energy);

  // FIXME: don't fire on the first inf
  if (std::abs(edrift) > delta_e_tol){
    std::cerr << "WARNING, energy drift exceeds tolerance!" << std::endl;
  }
  
  if (masses.has_inf() || arma::any(masses==0)){
    throw std::logic_error("Cannot do surface hopping with massless atoms or momentum sinks!");
  }

  // 3N vector of masses
  arma::vec m (3 * (NQM() + NMM()), arma::fill::zeros);
  for(arma::uword i = 0 ; i < m.n_elem ; i++){
    // no reason to do this without a bounds [] check!
    m(i) =  masses(i/3);
  }
    
  bool update = hopping; // copy so hopping can be reset
  double deltaE = 0;
  if (hopping){
    std::cerr <<  "[FSSH] Attempting hop: " << active_state + min_state
              << "->" << target_state + min_state << std::endl;
    
    deltaE = hop_and_scale(total_gradient, velocities, m);
    
    if (decoherence && deltaE > 0){
      decoherence->hopped(c, active_state);
    }
  }

  if (decoherence){
    decoherence->decohere(c, U, active_state, velocities.as_col(), m);
  }
  
  // update indicates we need to copy velocity and gradient back to GMX
  return update;
}
