#include "decoherence_afssh.hpp"

double hypot(std::complex<double> a, std::complex<double> b);

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

  dF_last.set_size(shstates);
  for (auto &v: dF_last) {v.set_size(3*(nqm + nmm)); v.zeros();}
  
  V.set_size(shstates);
}

/* Main entrance to decoherence routine*/
bool AFSSH::decohere(Electronic &c, const arma::mat U, const size_t active_state, const arma::vec v, const arma::vec m){
  bool collapsed = false;

  /* Jain step 7 */
  evolve_moments(c, U, active_state, m);

  /* Jain step 8 */

  // compute decoherence/reset rates
  arma::vec invtau_d(nstates); // collapse rate
  arma::vec invtau_r(nstates); // reset rate

  build_rates(invtau_d, invtau_r, U, active_state, v);

  // Collapse & reset as necessary
  {
    const arma::uword a = active_state;
    const double eta = arma::randu();
    for (arma::uword j = 0; j < nstates; j++){
      // Jain 8a
      if (eta < dtc * invtau_d(j)){
        c[a] = c(a) / std::abs(c(a)) * hypot(c(a), c(j));
        c[j] = 0;

        reset_moments(j);
        collapsed = true;

        std::cerr << "[AFSSH] Decoherence event: " << j << "->" << a << "."  << std::endl;
      }

      if (eta < dtc * invtau_r(j)){
        reset_moments(j);
        std::cerr << "[AFSSH] Reset " << j << " moments" << std::endl;
      }
    }
  }
  
  return collapsed;
}

void AFSSH::build_rates(arma::vec &invtau_d, arma::vec &invtau_r, const arma::mat U, const size_t a, const arma::vec v){
  const arma::mat T = arma::real(arma::logmat(U)) / dtc;
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
    double itau_fluct =  2 * arma::norm( T(a,n) * (V(a) - V(n)) * v.t() * (dR(n) - dR(a)) ) / arma::as_scalar(v.t() * v);
    invtau_d(n) = -1.0 * invtau_r(n) - itau_fluct;
    /*
      FIXME: Since MM forces are simply included in moment evolution,
      we will want to monitor the relative size of itau_fluct with MM
      system size.
    */
  }
}

/* Jain step 7 */
void AFSSH::evolve_moments(const Electronic &c, const arma::mat U, const size_t a, const arma::vec m){
  static bool first_call = true;
  
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
    props.emplace(QMProperty::energies, {nstates}, &energy);
    (*qm)->get_properties(props);

    // trim the energies that are not needed
    V = energy.tail(nstates);
    // V used in build_rates()
  }

  // gqm,gmm -> dF, a vector of length 3 (nqm+nmm) 
  for(arma::uword j = 0; j < nstates; j++){
    dF(j) = arma::vectorise(arma::join_horiz(gqm.slice(j) - gqm.slice(a),
                                             gmm.slice(j) - gmm.slice(a)));
    dF(j) *= -1.0;  // R.B.: we compute gradients, but need forces
  }

  // diabatic dF; eq 49
  for(arma::uword j = 0; j < nstates; j++){
    dF_(j).zeros();
    for(arma::uword k = 0; k < nstates; k++){
      dF_(j) += dF(k) * U(j,k) * U(j,k);
    }
  }

  if(first_call){
    save(dF_, c);
  }
  
  // evolve diabatic dR_ & dP_ moments; eqs 47, 48
  for(arma::uword j = 0; j < nstates; j++){
    dR_(j) += dtc * (dP_(j)  + 0.5*dtc*dF_(j)*rho(j)) / m ;
    dP_(j) += dtc * 0.5 * (dF_(j) + dF_last(j)) * rho(j);
  }  

  // Transform dR, dP moments: diabatic -> adiabatic; eqs 50, 51
  for(arma::uword j = 0; j < nstates; j++){
    dR(j).zeros();
    dP(j).zeros();
    for(arma::uword k = 0; k < nstates; k++){
      dR(j) += dR_(k) * U(k,j) * U(k,j);
      dP(j) += dP_(k) * U(k,j) * U(k,j);
    }
  }

  /*
    N.B.: The transformation taking dRj -> dRj - dRa
    is missing from Jain (2016). Amber confirmed that
    this is necessary in a call with DVCS on 2021.03.08
  */
  for(arma::uword j = 0; j < nstates; j++){
    dR(j) -= dR(a);
    dP(j) -= dP(a);
  }

  // store dF_ and rho for the next timestep 
  save(dF_, c);
  first_call = false;
}

/* Jain 2016 step 6; reset moments in event of hop */
void AFSSH::hopped(Electronic &c, size_t active_state){
  (void) c; (void) active_state;
  for (arma::uword j = 0; j < nstates; j++){
    reset_moments(j);
  }
}

/*
  std::hypot() for complex<double>

  returns sqrt(|a|^2 + |b|^2)

  I'm not sure what is the "best" sequence for the real() operations.
  The present choice---to call real() twice, once after each
  multiplication---was made because that's the first time the
  opperands mathematically lie on the real line. Waiting longer seemed
  to invite the accumulation of small numerical deviations.
*/
double hypot(std::complex<double> a, std::complex<double> b){
  return std::sqrt(std::real(std::conj(a)*a) + std::real(std::conj(b)*b) );
}

/* Jain 8b */
void AFSSH::reset_moments(arma::uword n){
  dR(n).zeros();  dR_(n).zeros();
  dP(n).zeros();  dP_(n).zeros(); 
}


// save dF and rho for use in next call.
void AFSSH::save(const arma::field<arma::vec> &dF_, const Electronic &c){
  dF_last = dF_;
  rho = arma::real(c.p().diag());
}
