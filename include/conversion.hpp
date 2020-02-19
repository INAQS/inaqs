#ifndef GIFS_CONVERSION_H
#define GIFS_CONVERSION_H

#include <algorithm>
#include <armadillo>


/* Direct conversion! */
class Conversion
{
public:
    Conversion(double in_en_from, double in_coords_to, double au_to_forces, double in_veloc_to) :
      en_from_au{in_en_from}, coords_to_au{in_coords_to}, veloc_to_au{in_veloc_to},
      au_to_forces {au_to_forces} //FIXME: MFSJM: please verify that this addition is as intended.
      {}

    /*
      FIXME: MFSJM: static, in this context, implies that there will
      only be one copy of this function shared between all instances
      of the Conversion class. As an attribute of a class member,
      static doesn't stay anything about linkage. Is this as intended?
    */
    static Conversion from_elementary(double mass, double length, double time);

    inline double energy_from_au(double en) { return en*en_from_au; }

    template<typename itr1, typename itr2>
    inline
    void 
    transform_coords_to_au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->coords_to_au;});
    };

    template<typename itr1, typename itr2>
    inline 
    void 
    transform_veloc_to_au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->veloc_to_au;});
    }

    template<typename itr1, typename itr2>
    inline 
    void 
    transform_au_to_forces(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->au_to_forces;});
    }
private:
    // energies
    double en_from_au;
    // coordinates
    double coords_to_au;
    // velocities
    double veloc_to_au;
    // forces
    double au_to_forces;
};

#endif
