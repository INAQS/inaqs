#include "decoherence_afssh.hpp"

// FIXME: should we include a way to specify the initial moments?
AFSSH::AFSSH(QMInterface ** const qm, const double dtc,
             const size_t min_state,
             const size_t shstates,
             const size_t nqm, const size_t nmm):
  Decoherence{qm, dtc, min_state, shstates, nqm, nmm} {
  // Set sizes for momments and potentials
  dR.set_size(shstates); dP.set_size(shstates); dF.set_size(shstates);
  for (auto &v: dR) {v.set_size(3*(nqm + nmm));}
  for (auto &v: dP) {v.set_size(3*(nqm + nmm));}
  for (auto &v: dF) {v.set_size(3*(nqm + nmm));}

  dR_.set_size(shstates); dP_.set_size(shstates); dF_.set_size(shstates);
  for (auto &v: dR_) {v.set_size(3*(nqm + nmm));}
  for (auto &v: dP_) {v.set_size(3*(nqm + nmm));}
  for (auto &v: dF_) {v.set_size(3*(nqm + nmm));}

  V.set_size(shstates);
}

/* Main entrance to decoherence routine*/
bool AFSSH::decohere(Electronic &c, const arma::mat U, const size_t active_state, const arma::vec v, const arma::vec m){
  std::cout << "AFSSH::decohere() on state " << active_state << std::endl;

  bool collapsed = false;

  /* Jain step 7 */
  evolve_moments(c, U, active_state, m);

  /* Jain step 8 */

  // compute decoherence/reset rates
  arma::vec invtau_d(nstates); // collapse rate
  arma::vec invtau_r(nstates); // reset rate

  build_invtau(invtau_d, invtau_r, U, active_state, v);

  {
    // N.B. we do this for states; multiple collapses update ca appropriately
    const double eta = arma::randu();
    for (arma::uword j = 0; j < nstates; j++){
      if (eta < dtc * invtau_d(j)){
        collapse(c, active_state, j);
        // FIXME: this should take the state as a parameter
        reset_moments();
        collapsed = true;
      }

      if (eta < dtc * invtau_r(j)){
        reset_moments();
      }
    }
  }
  
  return collapsed;
}

void AFSSH::build_invtau(arma::vec &invtau_d, arma::vec &invtau_r, const arma::mat U, const size_t a, const arma::vec v){
  const arma::mat T = arma::real(arma::logmat(U)) / dtc;
  // FIXME: (zeyu?) should invtau_d contain (-invtau_r) ?
  
  /* 
     Below, see an equivalent expression for tau that uses Ndim x
     Nstates matricies for moments. I don't think its as transparent
     as the for loop; for small numbers of states; I can't imagine the
     difference will matter.
  */
    
  // invtau_r = arma::sum(dF % (dR.each_col() - dR.col(a)))
    
  for (arma::uword n = 0; n < nstates; n++){
    // jain eq. 19
    invtau_r(n) = -1.0 * arma::as_scalar(dF(n).t() * (dR(n) - dR(a))) / 2;
    
    // jain eq. 54
    invtau_d(n) = -1.0 * invtau_r(n) -
      2 * arma::norm( T(a,n) * (V(a) - V(n)) * v.t() * (dR(n) - dR(a)) ) / arma::as_scalar(v.t() * v);
  }
}

/* Jain step 7 */
void AFSSH::evolve_moments(const Electronic &c, const arma::mat U, const size_t a, const arma::vec m){
  (void) U;
  (void) m;


  arma::cube gqm(3, nqm, nstates);
  arma::cube gmm(3, nmm, nstates);

  // This block fills gqm, gmm, and V with nstates objects
  {
    /*
      Two notes on this secion: 1) For all calls to
      (*qm)->get_properties(), we need to add min_state as in
      fssh.cpp. 2) There's no need to attempt to get the active state
      gradient first; qm_qchem will recalculate it, but will skip over
      scfman/setman so it'll be fast.
    */
    
    arma::uvec states(nstates);
    arma::vec energy(nstates + min_state);
  
    for (arma::uword i = 0; i < nstates; i++){
      states(i) = min_state + i;
    }
  
    PropMap props{};
  
    props.emplace(QMProperty::qmgradient_multi, states, &gqm);
    props.emplace(QMProperty::mmgradient_multi, states, &gmm);
    // R.B.: with any idx, energies will get all (excited_states + 1) states
    // FIXME: should make energies behavior consistent with gradients
    props.emplace(QMProperty::energies, {nstates}, &energy);
    (*qm)->get_properties(props);

    // trim the energies that are not needed
    V = energy.tail(nstates);
  }

  // Set dF
  for(arma::uword j = 0; j < nstates; j++){
    dF(j) = arma::vectorise(arma::join_horiz(gqm.slice(j) - gqm.slice(a),
                                             gmm.slice(j) - gmm.slice(a)));
  }

  // diabatic dF; eq 49
  for(arma::uword j = 0; j < nstates; j++){
    dF_(j).zeros();
    for(arma::uword k = 0; k < nstates; k++){
      dF_(j) += dF(k) * U(j,k) * U(j,k);
    }
  }


  /* 
     FIXME: Need to think carefully about when this is called and how
     to make sure we have the correct forces (at t0 and (t0 + dtc));
     it may be that we need to have mutiple hooks and conform to
     Gromacs' verlet-as-leapfrog implementation.
  */

  
  
  // evolve diabatic moments; eqs 47, 48
  for(arma::uword j = 0; j < nstates; j++){
    dR_(j) += dtc * (dP_(j)  + 0.5*dtc*dF_(j)*std::real(c.p(j,j))) / m ;
  }

  for (arma::uword j = 0; j < nstates; j++){
    // FIXME: need to do this in two passes with 0.5 F(t0) + 0.5 F(t0 + dtc);
    dP_(j) += dtc * dF_(j) * std::real(c.p(j,j));
  }
  

  // diabatic -> adiabatic R and P moments; eqs 50, 51
  for(arma::uword j = 0; j < nstates; j++){
    dR(j).zeros();
    dP(j).zeros();
    for(arma::uword k = 0; k < nstates; k++){
      dR(j) += dR_(k) * U(k,j) * U(k,j);
      dP(j) += dP_(k) * U(k,j) * U(k,j);
    }
  }

  // FIXME: (joe?) How to deal with mm region/forces in moment evolution?
  // FIXME: (joe?) From which classical time-step should we use the density matrix?  
}

/* Jain 2016 step 6; reset moments in event of hop */
void AFSSH::hopped(Electronic &c, size_t active_state){
  (void) c; (void) active_state;
  reset_moments();
}

/*
  std::hypot() for complex<double>

  returns sqrt(|a|^2 + |b|^2)

  I'm not sure what is the "best" sequence for the real() operations.
  The choice to call real() twice, once after each multiplication was
  made because that's the first time the opperands mathematically lie
  on the real line. Waiting longer seemed to invite the accumulation
  of small numerical deviations.
*/
double hypot(std::complex<double> a, std::complex<double> b){
  return std::sqrt(std::real(std::conj(a)*a) + std::real(std::conj(b)*b) );
}

/* Jain 8a */
void AFSSH::collapse(Electronic &c, size_t active_state, size_t collapse_state){
  const arma::uword a = active_state;
  const arma::uword n = collapse_state;

  std::complex<double> ca_new = c(a) / std::abs(c(a)) * hypot(c(a), c(n));
  c[a] = ca_new;
  c[n] = 0;
}

/* Jain 8b */
void AFSSH::reset_moments(void){
  for (auto &v: dR) {v.zeros();}  for (auto &v: dR_) {v.zeros();}
  for (auto &v: dP) {v.zeros();}  for (auto &v: dP_) {v.zeros();}
}
