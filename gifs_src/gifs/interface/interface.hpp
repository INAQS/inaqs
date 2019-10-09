#ifndef GIFS_SH_INTERFACE_CORE_H
#define GIFS_SH_INTERFACE_CORE_H

#include <string>
#include "gifs/molecule/system.hpp"

namespace gifs 
{
/* Virtual base class, defining all operations that 
 * the qm interface can perform
 */
class Interface
{
public:
    explicit Interface() {};
    virtual ~Interface() {};

    // compute functions
    // mechanical embedding:
    virtual unit get_energy(System& system) = 0;
    virtual unit get_eandg(System& system, std::vector<unit>& grad) = 0;
};


} // end namespace gifs

#endif // GIFS_SH_INTERFACE_CORE_H
