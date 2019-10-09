#ifndef GIFS_SH_GIFS_CORE_H
#define GIFS_SH_GIFS_CORE_H

#include <string>
#include <memory>

#include "gifs/sh/sh.hpp"
#include "gifs/es/es.hpp"

namespace gifs 
{
/* Virtual base class, defining all operations that 
 * should be called from the MD software
 *
 * P
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
    virtual unit get_energy(std::vector<double>& crd) = 0;
    virtual unit get_eandg(std::vector<double>& crd, std::vector<unit>& grad) = 0;
private:
    std::string _configfile{};
    // Implementation
    std::unique_prt<SurfaceHopping> _sh{nullptr};
    std::unique_prt<ElectronicStructure> _es{nullptr};
};


} // end namespace gifs

#endif // GIFS_SH_GIFS_CORE_H
