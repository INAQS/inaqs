#include "properties.hpp"
#include "constants.hpp"
#include "bomd.hpp"
#include "qm_qchem.hpp"
//
#include <armadillo>
//
QMInterface*
select_interface(ConfigBlockReader& reader,
                 FileHandle& fh,
                 const arma::uvec& qmids, 
	             arma::mat& qm_crd, 
	             arma::mat& mm_crd, 
	             arma::vec& mm_chg)
{
    std::string qmcode{};
    int mult, chg, nstates;
    reader.get_data("qm_code", qmcode);
    reader.get_data("mult", mult);
    reader.get_data("chg", chg);
    reader.get_data("nstates", nstates);

    if (qmcode == "qchem") {
        return new QM_QChem(fh, qmids, qm_crd, mm_crd, mm_chg, chg, mult, nstates);
    } else {
        throw "qm interface  not implemented!";
    }
}


BOMD::BOMD(arma::mat& qm_grd,
           arma::mat& mm_grd) :
    qm_grd{qm_grd}, mm_grd{mm_grd}, energy(1)
{};

void
BOMD::setup(FileHandle& fh,
            const arma::uvec& atomicnumbers,
            arma::mat& qm_crd, 
            arma::mat& mm_crd, 
            arma::vec& mm_chg
        ) {
  auto reader = setup_reader();
  add_necessary_keys(reader);
  reader.parse(fh);
  qm = select_interface(reader, fh, atomicnumbers, qm_crd, mm_crd, mm_chg);
  get_reader_data(reader);
};


void
BOMD::get_reader_data(ConfigBlockReader& reader) {
    const auto& reader2 = reader;
    (void) reader2;
};


void
BOMD::add_necessary_keys(ConfigBlockReader& reader) 
{
    reader.add_entry("qm_code", "qchem");
    reader.add_entry("chg", 0);
    reader.add_entry("mult", 1);
    reader.add_entry("nstates", 1);
};

ConfigBlockReader
BOMD::setup_reader() {
    //using types = ConfigBlockReader::types;
    ConfigBlockReader reader{"bomd"};
    //
    reader.add_entry("active_state", 0);
    return reader;
};

double 
BOMD::update_gradient()
{
    qm->update();

    PropMap props{};
    props.emplace(QMProperty::qmgradient, &qm_grd);
    props.emplace(QMProperty::mmgradient, &mm_grd);
    props.emplace(QMProperty::energies, &energy);
    //
    qm->get_properties(props);
    //
    return energy[0];
};


bool BOMD::rescale_velocities(arma::mat &velocities, arma::vec &masses, arma::mat &total_gradient, double total_energy) {
  (void) total_gradient;
  (void) masses;
  (void) velocities;
  edrift = (total_energy - elast)/elast;
  elast = total_energy;
  std::cout << "Total Energy: " << elast << ", " << edrift << std::endl;
  return false;
};
