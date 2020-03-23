#include <algorithm>
#include <armadillo>

#include "conversion.hpp"

// you own it!
Conversion* Conversion::from_elementary(double mass, double length, double time) {
    const double internal_length = 5.29177210903e-11;  // m, Bohr
    const double internal_mass =  9.1093837015e-31;    // kg, au
    const double internal_time =  2.418884326509e-17;  // seconds, au
    const double internal_velocity =  internal_length/internal_time;  // m/s
    const double internal_energy = internal_mass * internal_length * internal_length /
                                    (internal_time * internal_time);
    const double internal_force = internal_energy/internal_length;
    //
    const double energy = mass * length * length/(time * time);
    const double force = energy/length;
    const double velocity = length/time;
    //
    return new Conversion(internal_energy/energy, length/internal_length, internal_force/force, velocity/internal_velocity);
};
