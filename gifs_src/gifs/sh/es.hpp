#ifndef GIFS_SH_ES_CORE_H
#define GIFS_SH_ES_CORE_H

#include <string>



namespace gifs 
{

/* Abstract base class containing everything 
 * one needs to perform the qm calculation 
 *
 * Interface already knows:
 *     * NAtoms
 *     * NStates
 *     * atomids (QM)
 */
class QMInput 
{
public:
    using unit = double;
    explicit QMInput(unit* in_crd) : crd{in_crd} {}
protected:
    unit* crd = nullptr;
    // settings if one computes gradient or not etc.
};

/* Adding electrostatic embedding support */
class QMMM_Input:
    public QMInput
{
public:
    using unit = QMInput::unit;
    explicit QMMM_Input(unit* in_crd, int in_nmm, 
                        unit* in_mmcrd, unit* in_mmchg) :
        QMInput(in_crd), nmm{in_nmm}, mmcrd{in_mmcrd}, mmchg{in_mmchg} 
        {}
protected:
    int nmm{};
    unit* mmcrd = nullptr;
    unit* mmchg = nullptr;
};
/* Abstract base class containing everything 
 * qm calculation should return
 */
class QMout
{
public:
    using unit = double;
    explicit QMout(unit* in_energies, unit* in_grad) :
        energies{in_energies}, grad{in_grad} {}
protected:
    unit* energies = nullptr;
    unit* grad = nullptr;
};

class QMMM_out:
    public QMout
{
public:
    explicit QMMMout(unit* in_energies, unit* in_grad, unit* in_gradmm) :
        QMout{in_energies, in_grad}, gradmm{in_gradmm} {}
protected:
    unit* gradmm = nullptr;
};

/* Virtual base class, defining all operations that the qm interface can perform */
class ElectronicStructure 
{
public:
    explicit ElectronicStructure(int in_NQM, int in_NStates, 
                                 std::vector<int>& in_atomids) : 
        NQM{in_NQM}, NStates{in_NStates}, atomids{in_atomids} {};
    virtual ~ElectronicStructure() {};
    // compute functions
    // mechanical embedding
    virtual void compute(QMInput& qmin, QMout& qmout) = 0;
    // electrostatic embedding
    virtual void compute(QMMM_Input& qmin, QMMM_out& qmout) = 0;
protected:
    int NQM{};
    int NStates{};
    std::vector<int> atomids{};
};


} // end namespace gifs

#endif // GIFS_SH_ES_CORE_H
