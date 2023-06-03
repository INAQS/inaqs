#include "bomd_electronic.hpp"
#include "electronic.hpp"
#include "constants.hpp"
#include "util.hpp"
#include <armadillo>
#include <complex>
#include <cmath>


ConfigBlockReader ElectronicBomd::setup_reader()
{
    ConfigBlockReader reader{"bomd-electronic"};
    reader.add_entry("backpropagate", 0);
    reader.add_entry("amplitude_file", "cs.hdf5");
    {
      std::vector<std::complex<double>> cs {};
      reader.add_entry("amplitudes", cs);
    }

    return reader;
}


void ElectronicBomd::get_reader_data(ConfigBlockReader& reader) {
  reader.get_data("amplitude_file", amplitude_file);

  /* added in BOMD::add_qm_keys() */
  reader.get_data("min_state", min_state);

  {
    int in_bool;
    reader.get_data("backpropagate", in_bool);
    backpropagate = in_bool;
  }

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

  {
    std::vector<std::complex<double>> cs_in {};
    reader.get_data("amplitudes", cs_in);
    if (cs_in.size() > 0){
      if (cs_in.size() != (arma::uword) nstates){
        throw std::runtime_error("Number of Amplitudes do not match number of dynamics surfaces!");
      }
      c = Electronic(arma::cx_vec(cs_in)); 
    }
    else{
      arma::cx_mat initial(nstates, nstates, arma::fill::eye);
      c = Electronic(initial);
    }

    if (!c.normed()){
      throw std::runtime_error("Amplitudes are not normed; check your input!");
    }
  }
}

// This is our primary hook into the Gromacs (or other) MD loop
double ElectronicBomd::update_gradient(void){
  qm->update();

  // get gradients, energies, and overlap
  {
    PropMap props{};
    props.emplace(QMProperty::qmgradient, {min_state + active_state}, &qm_grd);
    props.emplace(QMProperty::mmgradient, {min_state + active_state}, &mm_grd);
    props.emplace(QMProperty::energies, util::range(min_state, min_state + nstates), &energy);
    props.emplace(QMProperty::wfoverlap, &U);
    qm->get_properties(props);
  }

  // write amplitudes first so we include the initial conditions
  c.get().save(arma::hdf5_name(amplitude_file,
                               "/amps/" + std::to_string(call_idx()),
                               arma::hdf5_opts::replace));
  shared->saveh5(U, "overlapraw");
  shared->saveh5(c.get(), "amps");

  Electronic::phase_match(U, phases);
  shared->saveh5(U, "overlap");
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
    if (!backpropagate){ // usual hamiltonian
      c.advance_exact(V - I*T, dtq);
    }
    else{
      c.advance_exact(-1.0*(V + I*T), dtq); 
    }
  }
    
  return energy(active_state);
}
