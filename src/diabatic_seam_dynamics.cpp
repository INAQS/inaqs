#include "diabatic_seam_dynamics.hpp"
#include <armadillo>
#include <utility>

ConfigBlockReader DiabaticSeam::setup_reader(){
  using types = ConfigBlockReader::types;
  ConfigBlockReader reader{"diabatic_seam"};
  
  reader.add_entry("adiabatI", types::ULINT);
  reader.add_entry("adiabatJ", types::ULINT);
  reader.add_entry("which_diabat", -1);
  
  return reader;
}

void DiabaticSeam::get_reader_data(ConfigBlockReader& reader) {
  reader.get_data("which_diabat", which_diabat);

  if (which_diabat != -1){
    single_diabat = true;

    if (!(which_diabat == 0 || which_diabat == 1)){
      throw std::runtime_error("Must specify 'which_diabat' as {0,1} to indicate donor or acceptor at starting geometry!");
    }
  }

  reader.get_data("adiabatI", lower);
  reader.get_data("adiabatJ", upper);

  size_t excited_states;
  /* added in BOMD::add_qm_keys() */
  reader.get_data("excited_states", excited_states);
  // R.B. logical OR has lower precedence than comparison
  // https://en.cppreference.com/w/cpp/language/operator_precedence
  if (lower > excited_states || upper > excited_states){
    throw std::runtime_error("Diabats out of the range of computed excited states; try increasing excited_states!");
  }

  if (upper < lower){
    std::swap(upper, lower);
  }
  
  energy.set_size(2);
  H.set_size(2,2);
  gd_qm.set_size(3,NQM(),2);
  g_qm.set_size(3,NQM(),2);
}


double DiabaticSeam::update_gradient(void){
  qm->update(); // Alert the QM module we have new data

  double E = 0;
  
  if (single_diabat){
    PropMap props {};
    props.emplace(QMProperty::diabatic_H, {lower, upper}, &H);
    props.emplace(QMProperty::diabatic_gradients, {lower,upper}, &gd_qm);
    qm->get_properties(props);

    if (call_idx() == 1){ // figure out which diabat is "upper" on the first step
      H.print("Initial diabatic Hamiltonian");
      std::cerr << "[DSEAM]: requested " << (which_diabat==0 ? "donor":"acceptor") << " diabat; "
                << "selecting state " << which_diabat << "." << std::endl;
    }

    shared->saveh5(H, "diabatic_H");

    E = H(which_diabat, which_diabat);
    qm_grd = gd_qm.slice(which_diabat);
  }
  else{
    PropMap props {};
    props.emplace(QMProperty::qmgradient_multi, {lower, upper}, &g_qm);
    props.emplace(QMProperty::energies, {lower,upper}, &energy);
    qm->get_properties(props);

    // FIXME: should save the diabatic Hamiltonian here too.
    
    E = 0.5 * (energy(0) + energy(1));
    qm_grd = 0.5 * (g_qm.slice(0) + g_qm.slice(1));
  }
  
  return E;
}
