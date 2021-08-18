#include <iostream>
#include <math.h>
#include "properties.hpp"
#include "constants.hpp"
#include "bomd_rescale.hpp"
#include "qm_qchem.hpp"
//
#include <armadillo>
//


void
RescaleBomd::get_reader_data(ConfigBlockReader& reader) {
//    reader.get_data("additional_energy", additional_energy);
//    reader.get_data("step", step);
//    reader.get_data("dE", dE);
    (void) reader;
};


ConfigBlockReader
RescaleBomd::setup_reader() {
    //using types = ConfigBlockReader::types;
    std::cout << "Setup Rescale BOMD\n";
    ConfigBlockReader reader{"rescalebomd"};
    //
    reader.add_entry("active_state", 0);
    /*
    reader.add_entry("additional_energy", 1.0);
    reader.add_entry("dE", 0.01);
    reader.add_entry("step", 1);
    */
    return reader;
};

double 
RescaleBomd::update_gradient()
{
    qm->update();

    /*
    PropMap props{};
    props.emplace(QMProperty::qmgradient, &qm_grd);
    props.emplace(QMProperty::mmgradient, &mm_grd);
    props.emplace(QMProperty::energies, &energy);
    //

    qm->get_properties(props);
    */
    //
    qm_grd.fill(0.0);
    mm_grd.fill(0.0);
    //
    double ene = 0.0; //energy[0];
    if (additional_energy > 0.0) {
       ene += additional_energy;
    }
    total_energy = ene;
    return ene;
};

double
kinetic_energy(arma::mat &velocities, arma::vec &masses) {
    double ekin = 0.0;
    for (arma::uword i=0; i<masses.n_elem; ++i) {
        for (arma::uword j=0; j<3; ++j) {
            ekin += 0.5 * velocities(j, i) * velocities(j, i)  * masses[i];
        }
    }
    return ekin;
}

bool RescaleBomd::rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double e_drift) {
  static int ncall = 1;
  bool retval = false;

  (void) total_gradient;
  (void) e_drift;

  const double ekin = kinetic_energy(velocities, masses); 
  const double factor = 1.0 + dE/ekin;
  //
  if ((ncall % step) == 0) {
    std::cout << "Rescale kinetic energy \n";
    std::cout << "factor: " << factor << "\n";
    if (additional_energy > dE) {
        std::cout << "we actually scale!\n";
        velocities *= sqrt(factor);
        additional_energy -= dE;
    }
    retval = true;
  }
  ncall += 1;
  std::cout << "Total energy " << total_energy << "au " << " ekin = " << ekin << "\n";
  return retval;
};
