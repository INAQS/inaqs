#include <algorithm>
#include <armadillo>

#include "conversion.hpp"

Conversion::Conversion(const double mass, const double length, const double time) {
    const double internal_length = 5.29177210903e-11;  // m, Bohr
    const double internal_mass =  9.1093837015e-31;    // kg, au
    const double internal_time =  2.418884326509e-17;  // seconds, au
    const double internal_velocity =  internal_length/internal_time;  // m/s
    const double internal_energy = internal_mass * internal_length * internal_length /
                                    (internal_time * internal_time);
    const double internal_gradient = internal_energy/internal_length;
    //
    const double energy = mass * length * length/(time * time);
    const double gradient = energy/length;
    const double velocity = length/time;
    //
    _energy_au2md = internal_energy/energy;
    _energy_md2au = energy/internal_energy;
    // coordinates
    _crd_au2md = internal_length/length;
    _crd_md2au = length/internal_length;
    // velocities
    _veloc_au2md = internal_velocity/velocity;
    _veloc_md2au = velocity/internal_velocity;
    // gradient
    _grd_au2md = internal_gradient/gradient;
    _grd_md2au = gradient/internal_gradient;
    // mass
    _mass_au2md = internal_mass/mass;
    _mass_md2au = mass/internal_mass; 
    //
};
