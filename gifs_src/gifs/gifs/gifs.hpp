#ifndef GIFS_SH_GIFS_CORE_H
#define GIFS_SH_GIFS_CORE_H

#include <string>
#include <memory>

namespace gifs 
{
/* Virtual base class, defining all operations that 
 * should be called from the MD software
 *
>>>>>>> 2e157ee07321291b7fe280d48e01e62c1363ef85
 */
class Gifs
{
public:
    using unit = double;
    // factory to generate SurfaceHopping code
    explicit Gifs(std::string fname);
    ~Gifs() {};
    // call surface hopping Routines
    void do_surface_hopping();
    // compute functions
    // mechanical embedding:
    unit get_energy(std::vector<unit>& crd);
    unit get_eandg(std::vector<unit>& crd, std::vector<unit>& grad);
private:
    std::string _configfile{};
    // Implementation
    std::unique_ptr<SurfaceHopping> _sh{nullptr};
};

} // end namespace gifs

#endif // GIFS_SH_GIFS_CORE_H
