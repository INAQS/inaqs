#include <vector>
#include "gifs.hpp"
#include "properties.hpp"
#include "qm.hpp"

static QMInterface* qm = nullptr;

void create_qm_interface(size_t nqm, const int* qm_atomids)
{
    if (qm == nullptr) {
        std::vector<int> qmids(qm_atomids, qm_atomids+nqm);
        qm = new QMInterface(nqm, qmids);
    };
}


float gifs_get_forces(const float* qm_crd, size_t nmm, const float* mm_crd, const float* mm_chg, float* f, float* fshift)
{
    size_t nqm = qm->get_nqm();

    qm->update(qm_crd, nmm, mm_crd, mm_chg);
    
    std::vector<double> g_qm(nqm*3), g_mm(3*nmm);
    std::vector<double> energy(1);

    PropMap props{};
    props.emplace(QMProperty::qmgradient, &g_qm);
    props.emplace(QMProperty::mmgradient, &g_mm);
    props.emplace(QMProperty::energies, &energy);
    //
    qm->get_properties(props);

    /* Unit conversion back qm->mm*/
    {
      int i = 0, j = 0;

      for (auto& fqm: g_qm) {
        j = i++;
        f[j] = HARTREE_BOHR2MD*fqm;
        fshift[j] = f[j];
      }

      for (auto& fmm: g_mm) {
        j = i++;
        f[j] = HARTREE_BOHR2MD*fmm;
        fshift[j] = f[j];
      }
    }

    /* Return in "MD" units of KJ/mole */
    return HARTREE2KJ * AVOGADRO * energy[0];

};
