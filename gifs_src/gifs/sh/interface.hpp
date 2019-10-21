#ifndef GIFS_SH_SURFACEHOPPING_CORE_H
#define GIFS_SH_SURFACEHOPPING_CORE_H

#include <string>
#include "gifs/molecule/system.hpp"

namespace gifs 
{
/* Virtual base class, defining all operations that 
 * the surface hopping code should provide
 */
class SurfaceHopping 
{
public:
    explicit SurfaceHopping() {};
    virtual ~SurfaceHopping() {};

    // compute functions
    // mechanical embedding:
    virtual unit get_energy(System& system) = 0;
    virtual unit get_eandg(System& system, std::vector<unit>& grad) = 0;

    virtual void do_surface_hopping()  = 0;
    // general 
    virtual double get_state_energy() = 0;
    // embedding dependent
    virtual double* get_state_gradient() = 0;
    int get_current_state() { return _istate; }
protected:
    int _istate{};
private:
};


} // end namespace gifs

#endif // GIFS_SH_SURFACEHOPPING_CORE_H
