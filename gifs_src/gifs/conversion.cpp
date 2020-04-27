#include <algorithm>
#include <armadillo>

#include "conversion.hpp"


Conversion* Conversion::from_elementary(const double mass, const double length, const double time) {
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
    return new Conversion(internal_energy/energy, internal_length/length, internal_gradient/gradient, internal_velocity/velocity, internal_mass/mass);
};
