#ifndef GIFS_SH_GIFS_CORE_H
#define GIFS_SH_GIFS_CORE_H

#include <string>
#include <memory>

namespace gifs 
{
/* Virtual base class, defining all operations that 
 * should be called from the MD software */
class Gifs
{
public:
    using unit = double;
    // factory to generate SurfaceHopping code
    explicit Gifs(std::string fname);
    ~Gifs() {};
    // call surface hopping Routines
    void do_surface_hopping();
    // rescaling (QM only? or also MM parts?)
    void rescale_velocities(unit* veloc, unit* mass);
    void rescale_velocities(int* iqm, unit* veloc, unit* mass);
    // compute functions
    // mechanical embedding:
    unit get_energy(unit* crd);
    unit get_eandg(unit* crd, unit* grad);
private:
    std::string _configfile{};
    // Implementation
    std::unique_ptr<SurfaceHopping> _sh{nullptr};
};

} // end namespace gifs

#endif // GIFS_SH_GIFS_CORE_H
