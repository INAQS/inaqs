#ifndef GIFS_SH_GIFS_CORE_H
#define GIFS_SH_GIFS_CORE_H

#include <string>
#include <memory>

namespace gifs 
{
/* Base class, defining all operations that 
 * should be called from the MD software/nuclear propagator
 */
class Gifs
{
public:
    explicit Gifs(std::string fname) : _configfile{fname} {};
    virtual ~Gifs() {};
    // surface hopping Routines
    virtual void do_surface_hopping()  = 0;
    // compute functions
    // mechanical embedding:
    inline double get_energy(std::vector<double>& crd) = 0;
    inline double get_eandg(std::vector<double>& crd, std::vector<unit>& grad) = 0;
private:
    std::string _configfile{};
};


} // end namespace gifs

#endif // GIFS_SH_GIFS_CORE_H
