#ifndef GIFS_SH_SH_CORE_H
#define GIFS_SH_SH_CORE_H

#include "gifs/sh/es.hpp"

namespace gifs 
{
/* Virtual base class, defining all operations that 
 * should be implemented in surface hopping software
 */
class SurfaceHopping
{
public:
    using unit = double;
    // setup specific SurfaceHopping code
    explicit SurfaceHopping(std::string input);
    ~SurfaceHopping() {};
    // virtual interface to call surface hopping Routines
    // especially electronic propagation etc.
    virtual void do_surface_hopping() = 0;
    // getter functions for the md code
    virtual unit get_energy(std::vector<unit>& crd) = 0;
    virtual unit get_eandg(std::vector<unit>& crd, std::vector<unit>& grad) = 0;
    // perform computations using electronic structure interface
    void compute(QMInput& qmin, QMout& qmout) {
        _es->compute(qmin, qmout);   
    };
    // electrostatic embedding
    void compute(QMMM_Input& qmin, QMMM_out& qmout) {
        _es->compute(qmin, qmout);
    };
protected:
    // Implementation of the electronic structure routines
    std::unique_ptr<ElectronicStructure> _es{nullptr};
};


} // end namespace gifs

#endif // GIFS_SH_SH_CORE_H
