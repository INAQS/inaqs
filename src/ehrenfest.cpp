#include "ehrenfest.hpp"
#include "electronic.hpp"
#include "constants.hpp"
#include "util.hpp"
#include <armadillo>
#include <complex>
#include <cmath>


ConfigBlockReader Ehrenfest::setup_reader()
{
    using types = ConfigBlockReader::types;
    ConfigBlockReader reader{"ehrenfest"};
    reader.add_entry("amplitude_file", "cs.dat");
    reader.add_entry("dtc", types::DOUBLE);

    {
      std::vector<std::complex<double>> cs {};
      reader.add_entry("amplitudes", cs);
    }

    return reader;
}


void Ehrenfest::get_reader_data(ConfigBlockReader& reader) {
  {
    double in_dtc;
    reader.get_data("dtc", in_dtc);  // in fs
    dtc = in_dtc * (1e-15 / AU2SI_TIME);  // fs -> a.u.
  }
  reader.get_data("amplitude_file", amplitude_file);

  /* added in BOMD::add_qm_keys() */
  reader.get_data("min_state", min_state);

  {
    size_t excited_states;
    /* added in BOMD::add_qm_keys() */
    reader.get_data("excited_states", excited_states);

    nstates = excited_states + 1 - min_state;
  }

  if (nstates < 2){
    throw std::logic_error("Doesn't make sense to propagate the electronic wavefunction on only 1 surface!");
  }

  if (!(min_state <= active_state && active_state <= min_state + nstates)){
    throw std::range_error("Active state not in the range of propagating states!");
  }

  active_state -= min_state;

  energy.set_size(nstates);

  U.set_size(nstates, nstates);
  T.set_size(nstates, nstates);
  V.set_size(nstates, nstates);
  phases.set_size(nstates);
  phases.ones();

  qm_grds.set_size(3, NQM(), nstates);

  {
    std::vector<std::complex<double>> cs_vec {};
    reader.get_data("amplitudes", cs_vec);
    if (cs_vec.size() > 0){
      if (cs_vec.size() != (arma::uword) nstates){
        throw std::runtime_error("Number of Amplitudes do not match number of hopping surfaces!");
      }
      double norm = 0;
      for (const auto& c: cs_vec){
        norm += std::norm(c);
      }
      std::cerr << "[Ehrenfest] amplitude norm = " << norm << std::endl;
      if (std::abs(1.0 -norm) > 1e-6){
        throw std::runtime_error("Amplitudes are not normed; check your input!");
      }
      arma::cx_vec cs_arma(nstates);
      for (int i = 0; i < nstates; i++){
        cs_arma(i) = cs_vec[i];
      }
      c = Electronic(cs_arma);
    }
    else{
      c.reset(nstates, 1, active_state);
    }
  }
}

// This is our primary hook into the Gromacs (or other) MD loop
double Ehrenfest::update_gradient(void){
  qm->update();
  mm_grds.set_size(3, NMM(), nstates);
  
  // get gradients, energies, and overlap
  {
    PropMap props{};
    auto states = util::range(min_state, min_state + nstates);
    props.emplace(QMProperty::qmgradient_multi, states, &qm_grds);
    props.emplace(QMProperty::mmgradient_multi, states, &mm_grds);
    props.emplace(QMProperty::energies, states, &energy);
    props.emplace(QMProperty::wfoverlap, &U);
    qm->get_properties(props);
  }
  
  // write populations/amplitudes first so we include the initial conditions
  {
    std::ofstream output(amplitude_file, std::ios_base::app);
    auto a = c();
    a.transform( [](std::complex<double> ci) { return std::norm(ci); } );
    arma::real(a).st().print(output);
    output.close();
  }
  saveh5(c.get(), "amps");  
  saveh5(U, "overlapraw");
  Electronic::phase_match(U, phases);
  saveh5(U, "overlap");
  T = util::logmat_unitary(U) / dtc;
  V = arma::diagmat(energy);

  // compute max time step for electronic propagation (Jain 2016 eqs. 20, 21)
  {
    double dtq_ =
      std::min(dtc/20, // make sure 20dtq < dtc
               std::min(dtc,
                        std::min(0.02 / T.max(),
                                 0.02 / arma::max( V.diag() - arma::mean(V.diag()) )
                                 )
                        ));

    dtq = dtc / std::round(dtc / dtq_);
  }

  const size_t n_steps = (size_t) dtc / dtq;
  const std::complex<double> I(0,1);

  // Propagate electronic coefficients
  for (size_t nt = 0; nt < n_steps; nt++){
    c.advance_exact(V - I*T, dtq);
  }

  double eh_energy = 0;
  mm_grd.zeros();
  qm_grd.zeros();
  for (int i = 0; i < nstates; i++){
    double ai = std::norm(c(i));
    eh_energy += ai * energy(i);
    qm_grd += ai * qm_grds.slice(i);
    mm_grd += ai * mm_grds.slice(i);
  }


  //coupling terms
  {
    arma::mat d(3, NQM()+NMM());
    // Upper triangle of dij
    for (arma::uword j=1; j < (arma::uword) nstates; j++){
      for (arma::uword i=0; i < j; i++){
        PropMap props{};
        props.emplace(QMProperty::nacvector, {min_state + i, min_state + j}, &d);
        d *= phases(i) * phases(j);
        qm->get_properties(props);
        const arma::mat & dqm = d.submat(arma::span::all, arma::span(0, NQM()-1));
        double cdv = 2 * (std::conj(c(i))*c(j)).real() * (energy(j) - energy(i));
        //std::cout << "[Ehrenfest] " << call_idx() << ": " <<
        //i << "," << j << ": " << cdv << std::endl;
        qm_grd += cdv * dqm;
        if (NMM() > 0){
          const arma::mat & dmm = d.submat(arma::span::all, arma::span(NQM(), NMM()-1));
          mm_grd += cdv * dmm;
        }
      }
    }
  }
  return eh_energy;
}
