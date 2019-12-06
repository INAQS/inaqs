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
    // 
    double getGradient(

            double* Grd
            
            
    );
    // rescaling (QM only? or also MM parts?)
    void rescale_velocities(unit* veloc, unit* mass, unit* gradient);
private:
    std::string _configfile{};
    // fixed
    int NQM;             // const
    int* atomids;        // NQM
    // 
    int* qm_idx;         // NQM, are they const throughout the simulation?
    // flexible
    int NMM;             // flexible
    double* crd_qm;      // NQM*3 
    double* crd_mm;      // NMM*3
    double* chg_mm;      // NMM
    // Implementation
    std::unique_ptr<SurfaceHopping> _sh{nullptr};
};

} // end namespace gifs

#endif // GIFS_SH_GIFS_CORE_H
