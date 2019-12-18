#ifndef __DYNAMICS_HPP
#define __DYNAMICS_HPP
#include <unordered_map>
#include <vector>
#include "properties.hpp"

/*
 * Basic strategie:
 *
 *    - Everything that is runtime changeable is handled via polymorphism
 *    - Factories are used to handle construction
 *
 * Main Classes:
 *    - BOMD:
 *      Handels actual requests
 *    - QMInterface:
 *      Interface to QM Software to perform actual computations
 *    - System
 *      Everything related to the qm/mm system,
 *      handels partitining of the forces etc.
 *    - Electronic:
 *      contains everything related to the wavefunction
 *      knows how to propagate and get back coefficients
 *      (probably also phasetracking etc.)
 *
 */
class BOMD
{
public:

    template<typename T>
    T get_gradient(T* qm_crd, T* mm_crd,
                   T* mm_chg, T* qm_gradient,
                   T* mm_gradient)
    {
        static_assert<std::is_floating_point<T>::value, "Error Msg">;
        ...
    }

    //
    template<typename T>
    inline
    T rescale_velocities(T* total_gradient, T* masses, T* velocities) { }

protected:
    QMInterface* QM;
    // fixed size
    std::vector<double> qm_grd;
    // flexible size
    std::vector<double> mm_grd;
};


class FSSH:
    public BOMD
{
    template<typename T>
    T get_gradient(T* qm_crd, T* mm_crd,
                   T* mm_chg, T* qm_gradient,
                   T* mm_gradient)
    {
        static_assert<std::is_floating_point<T>::value, "Error Msg">;
        ...
        // propagate electronic wavefunction
        QM-get_properties();
    }

    //
    template<typename T>
    inline
    T rescale_velocities(T* total_gradient, T* masses, T* velocities) 
    {
        static_assert<std::is_floating_point<T>::value, "Error Msg">;
        // get local qm representation of the properties
        sys->from_global_to_qm(total_gradient, qm_total_grad);
        ...
        //
        if (is_leapfrog) {
            leapfrog_backwards(qm_total_grad, qm_masses, qm_velocities);
        }

        double ekin = compute_ekin(qm_masses.size, qm_masses, qm_velocities);

        /*
         * ATTENTION:
         *         - This is currently wrong.....
         *         - should be in get_gradient....
         *         - actually should be somewhere in between....
         *         - if we want to handel trivial hops, it gets somehow complicated...
         *
         * SOLUTION:
         *         Default: do not handle trivial hops
         *         Else: ensure no constrains on the QM region!
         *
         *
         */
        if (try_hopping()) {
            // update qm_grd and mm_grd and potentially compute stuff
            QM->get_properties();
            //
            rescale_vector(scale, qm_velocities);
            //
            if (is_leapfrog) {
                leapfrog_forward(qm_total_grad, qm_masses, qm_velocities);
            }
            //
            /* copy data back */
            sys->from_qm_to_global(qm_total_grad, total_gradient);
            ...
            // this makes no sense at the current position!
            if (is_trivial_hop) {
                // update gradient
                update_gradient(qm_total_grad);
                // handle LA
                sys->from_qm_to_global(qm_total_grad, total_gradient);
            } else {
                ...
            }
        }
        // do the decoherence correction
        decoherence();
    }


protected:
    //
    bool try_hopping() {
        // everything related to the surface hopping algorithm
        bool hopped = false;
        for (int isubstep=0; isubstep < nsubsteps; ++isubstep) {
            // handle somehow linear interpolation for H
            Psi->propagate();
            //
            if (!hopped) {
                coefs = Psi->get_coefficients();
                compute_hopping(coefs)
            }
       }

    }

    void decoherence() {
        //
    }

    //
    void update_gradient() {
    }
    //
    Electronic* Psi;
    System* sys;
    //
    std::vector<double> energies; // NSTATES
    std::vector<double> overlaps; // NSTATES*NSTATES
    // only necessary if rescale via nacs
    std::vector<double> nac;  // 3*NQM
    //
    bool is_leapfrog{false};
    // info about bounderies!
};

/*
 * Contains all info regarding the electronic wavepacket!
 * can give back coefficents (needed for surface hopping)
 * and handels electronic propagation
 *
 */
class Electronic
{
public:
    // propagate electronic wavefunction
    void propagate();
    void get_coefficients(int wfid);
    // get first entry
    void get_coefficients();
    //
private:
};


/*
 * Class to handle all possible linkatom
 * schemes, is very much implementation dependent
 * and therefore the gifs library user should be
 * able to define there own!
 */
class System
{
public:
    // get from a global system object, the local representation
    void from_global_to_qm(total, qm);
    //
    void from_qm_to_global(qm, total);
    // add/subtract qm/mm forces
    void subtract_forces(qm, mm, total);
    void add_forces(qm, mm, total);
    // Should also handle the non link atom case, where only qm_idx are used!
protected:
    int NLink;                        // const
    std::vector<int> qm_idx;         // NQM, are they const throughout the simulation?
    // QM-MM
    linkatoms array, (idx, QM, MM)      //
};

/*
 * Linkatoms are treated outside,
 * the qminterface class does not know how to handle them
 * and not even that they exist!
 * same with MM atoms, and fake MM atoms (rcD etc.)
 *
 */
class QMInterface
{
public:
    QMInterface(int nqm, std::vector<int> qmid};
    void get_properties(PropMap &props);

    void update(std::vector<double> &crdqm,
		std::vector<double> &crdmm,
		std::vector<double> &chgmm);  

protected:
    int NQM;             // const, actually NQM+NLink
    //  fixed size
    std::vector<int> atomids;        // NQM
    std::vector<double> crd_qm;      // NQM*3
    // flexible
    int NMM;
    std::vector<double> crd_mm;      // NMM*3
    std::vector<double> chg_mm;      // NMM

};

#endif
