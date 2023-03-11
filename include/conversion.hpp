#ifndef GIFS_CONVERSION_H
#define GIFS_CONVERSION_H

#include <algorithm>
#include <armadillo>


/* Direct conversion! */
class Conversion
{
public:
    // 
    explicit Conversion(double mass, double length, double time);
    Conversion(const Conversion& rhs) = default;
    Conversion(Conversion&& rhs) = default;
    Conversion& operator=(const Conversion& rhs) = default;
    Conversion& operator=(Conversion&& rhs) = default;

    //
    template<typename itr1, typename itr2>
    inline
    void 
    transform_crd_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double value){ return this->crd_md2au(value); });
    };
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_veloc_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double value){ return this->veloc_md2au(value); });
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_veloc_au2md(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double value){ return this->veloc_au2md(value); });
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_gradient_au2md(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double value){ return this->grd_au2md(value); });
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_gradient_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double value){ return this->grd_md2au(value); });
    }
    //
    template<typename itr1, typename itr2>
    inline 
    void 
    transform_masses_md2au(itr1 in_begin, itr1 in_end, itr2 result_begin) {
        std::transform(in_begin, in_end, result_begin, [this](double value){ return this->mass_md2au(value); });
    }
    //
    inline double energy_au2md(double value) const noexcept {return _energy_au2md*value;}
    inline double energy_md2au(double value) const noexcept {return _energy_md2au*value;}
    //
    inline double crd_au2md(double value) const noexcept {return _crd_au2md*value;}
    inline double crd_md2au(double value) const noexcept {return _crd_md2au*value;}
    //
    inline double veloc_au2md(double value) const noexcept {return _veloc_au2md*value;}
    inline double veloc_md2au(double value) const noexcept {return _veloc_md2au*value;}
    //
    inline double grd_au2md(double value) const noexcept {return _grd_au2md*value;}
    inline double grd_md2au(double value) const noexcept {return _grd_md2au*value;}
    //
    inline double mass_au2md(double value) const noexcept {return _mass_au2md*value;}
    inline double mass_md2au(double value) const noexcept {return _mass_md2au*value;}
      //
    inline double time_au2md(double value) const noexcept {return _time_au2md*value;}
    inline double time_md2au(double value) const noexcept {return _time_md2au*value;}

private:
    // energies
    double _energy_au2md;
    double _energy_md2au;
    // coordinates
    double _crd_au2md;
    double _crd_md2au;
    // velocities
    double _veloc_au2md;
    double _veloc_md2au;
    // gradient
    double _grd_au2md;
    double _grd_md2au;
    // mass
    double _mass_au2md;
    double _mass_md2au;
    // time
    double _time_au2md;
    double _time_md2au;
};

#endif
