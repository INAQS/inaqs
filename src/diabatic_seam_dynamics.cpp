#include "diabatic_seam_dynamics.hpp"
#include <armadillo>

ConfigBlockReader DiabaticSeam::setup_reader(){
  using types = ConfigBlockReader::types;
  ConfigBlockReader reader{"diabatic_seam"};
  
  reader.add_entry("adiabatI", types::ULINT);
  reader.add_entry("adiabatJ", types::ULINT);
  reader.add_entry("single_diabat", 0);
  
  return reader;
}

void DiabaticSeam::get_reader_data(ConfigBlockReader& reader) {
  {
    int bool_in;
    
    reader.get_data("single_diabat", bool_in);
    single_diabat = bool_in;
  }

  {
    size_t I, J;
    reader.get_data("adiabatI", I);
    reader.get_data("adiabatJ", J);

    size_t excited_states;
    /* added in BOMD::add_qm_keys() */
    reader.get_data("excited_states", excited_states);
    // R.B. logical OR has lower precedence than comparison
    // https://en.cppreference.com/w/cpp/language/operator_precedence
    if (I > excited_states || J > excited_states){
      throw std::runtime_error("Diabats out of the range of computed excited states; try increasing excited_states!");
    }

    //if (I == J){
    // throw std::runtime_error("Must specify different adiabats to diabatize!");
    //}

    if (I > J){
      upper = I;
      lower = J;
    }
    else{
      upper = J;
      lower = I;
    }
  }
  
  energy.set_size(2);
  diabatic_rot_mat.set_size(2,2);
  gd_qm.set_size(3,NQM(),2);
  g_qm.set_size(3,NQM(),2);
}


double DiabaticSeam::update_gradient(void){
  /*
    FIXME: Allow running on single diabat
  */
  
  qm->update();
  PropMap props {};
  props.emplace(QMProperty::qmgradient_multi, {lower, upper}, &g_qm);
  props.emplace(QMProperty::energies, {lower,upper}, &energy);
  qm->get_properties(props);

  double E = 0.5 * (energy(0) + energy(1));
  qm_grd = 0.5 * (g_qm.slice(0) + g_qm.slice(1));

  return E;
  
  
  // //FIXME: probably want to expose dipoles to this class
  // PropMap props{};

  // props.emplace(QMProperty::qmgradient, {upper}, &g_qm.slice(1));  
  // props.emplace(QMProperty::diabatic_gradients, {lower, upper}, &gd_qm);  
  // props.emplace(QMProperty::diabatic_rot_mat, {lower, upper}, &diabatic_rot_mat);

  // arma::vec absolute_energies(3);
  // props.emplace(QMProperty::energies, {0,lower,upper}, &absolute_energies);
  // qm->get_properties(props);

  // // need energies to be excitations
  // energy = absolute_energies.subvec(1,2) - absolute_energies(0);
  
  // diabatic_energy = (diabatic_rot_mat.t() * arma::diagmat(energy) * diabatic_rot_mat).eval().diag();
  
  // saveh5(diabatic_energy, "diabats");
  // saveh5(energy, "adiabats");
  // saveh5(diabatic_rot_mat, "rotation");
  
  // std::cerr << "[DSD] " << call_idx() << ": gap = " << (diabatic_energy(0) - diabatic_energy(1)) << std::endl;  

  // if (single_diabat){
  //   qm_grd = gd_qm.slice(0);
  //   return diabatic_energy(0);
  // }
  
  // if (alpha != 0){
  //   return build_diabatic_forces_restrained();
  // }
  // else{
  //   return build_diabatic_forces_projected();
  // }
  qm->update();
}


// FIXME: this projection scheme cannot conserve energy; figure out
// why analytically and rectify if possible. (Might require projecting
// out a component of the velocities too in rescale_velocities.)
double DiabaticSeam::build_diabatic_forces_projected(void){
  double E = diabatic_energy(0);

  arma::vec v = gd_qm.slice(0).as_col() - gd_qm.slice(1).as_col();
  v = arma::normalise(v);

  
  arma::vec projection = v * arma::as_scalar(v.t() * gd_qm.slice(0).as_col());  
  qm_grd = gd_qm.slice(0) - arma::reshape(projection, arma::size(qm_grd));
  
  return E;
}


double DiabaticSeam::build_diabatic_forces_restrained(void){
  double delta = (diabatic_energy(0) - diabatic_energy(1)); 
  double restraint = 0.5 * alpha * delta*delta;
  std::cerr << "[DSD] " << call_idx() << ": restraint = " << restraint << std::endl;  
  double E = energy(1) + restraint; // take the upper *adiabat* + restraint
  
  qm_grd = g_qm.slice(1) + alpha*delta*(gd_qm.slice(0) - gd_qm.slice(1));

  return E;

  /*
    double delta = (diabatic_energy(0) - diabatic_energy(1)); 
    double restraint = 0.5 * alpha * delta*delta;
    std::cerr << "[DSD] " << call_idx() << ": restraint = " << restraint << std::endl;  
    double E = 0.5 * (diabatic_energy(0) + diabatic_energy(1)) + restraint;

    std::cerr << "[DSD] " << call_idx() << ": |F| = " << arma::norm((gd_qm.slice(0) + gd_qm.slice(1)).as_col()) <<
    ", |R| = " << arma::norm((alpha*delta*(gd_qm.slice(0) - gd_qm.slice(1))).as_col()) << std::endl;  
  
    qm_grd = 0.5 * (gd_qm.slice(0) + gd_qm.slice(1)) + alpha*delta*(gd_qm.slice(0) - gd_qm.slice(1));
  
    return E;
  */
}
