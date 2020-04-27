#ifndef GIFS_CONVERSION_H
#define GIFS_CONVERSION_H

#include <algorithm>
#include <armadillo>


/* Direct conversion! */
class Conversion
{
public:
    // public constructor
    static Conversion* from_elementary(double mass, double length, double time);
    //
    template<typename itr1, typename itr2>
    inline
    void 
    transform_coords_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->crd_md2au;});
    };
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_veloc_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->veloc_md2au;});
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_veloc_au2md(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->veloc_au2md;});
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_gradient_au2md(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->grd_au2md;});
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_gradient_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->grd_md2au;});
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_masses_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double val) -> double {return val*this->mass_md2au;});
    }
    // energies
    const double energy_au2md;
    const double energy_md2au;
    // coordinates
    const double crd_au2md;
    const double crd_md2au;
    // velocities
    const double veloc_au2md;
    const double veloc_md2au;
    // gradient
    const double grd_au2md;
    const double grd_md2au;
    // mass
    const double mass_au2md;
    const double mass_md2au;
private:
    explicit Conversion(
            double in_energy_md2au, 
            double in_crd_md2au, 
            double in_grd_md2au, 
            double in_veloc_md2au, 
            double in_mass_md2au) :
    energy_au2md{1.0/in_energy_md2au},
    energy_md2au{in_energy_md2au},
    // coordinates
    crd_au2md{1.0/in_crd_md2au},
    crd_md2au{in_crd_md2au},
    // velocities
    veloc_au2md{1.0/in_veloc_md2au},
    veloc_md2au{in_veloc_md2au},
    // gradient
    grd_au2md{1.0/in_grd_md2au},
    grd_md2au{in_grd_md2au},
    // mass
    mass_au2md{1.0/in_mass_md2au},
    mass_md2au{in_mass_md2au}
    {}
};

#endif
