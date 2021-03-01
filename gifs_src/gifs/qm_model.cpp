#include <ctime>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <algorithm>
#include "qm_model.hpp"
#include "properties.hpp"
#include <armadillo>
#include <unordered_map>

ConfigBlockReader QM_Model::model_reader() {
  ConfigBlockReader reader{"qmmodel"};
  reader.add_entry("model", "avoided-crossing");
  return reader;
}

QM_Model::QM_Model(FileHandle& fh, 
                   const arma::uvec& in_qmids, 
		   arma::mat& in_qm_crd, 
		   arma::mat& in_mm_crd, 
		   arma::vec& in_mm_chg, 
		   const int charge, 
		   const int mult,
		   const int excited_states,
                   const int min_state):
  QMInterface(in_qmids, in_qm_crd, in_mm_crd, in_mm_chg, charge, mult, excited_states, min_state)
{

  if (NQM > 1){
    throw std::logic_error("Model systems only work on the x component of a singele atom!");
  }

  if (NMM != 0){
    throw std::logic_error("Model systems don't support QM/MM!");
  }
  
  ConfigBlockReader reader = model_reader();
  reader.parse(fh);

  {
    std::string model_in;
    reader.get_data("model", model_in);

    if (model_in == "avoided-crossing"){
      model= new AvoidedCrossing();
    }
    else{
      throw std::runtime_error("Unknown Model: `" + model_in + "`!");
    }
  }
}

void QM_Model::update(void){
  QMInterface::update();
  // always take the x component of the first atom
  model->update(crd_qm(0));
}

void QM_Model::get_properties(PropMap &props){
  for (QMProperty p: props.keys()){
    switch(p){
      
    case QMProperty::wfoverlap:{
      arma::mat &U = props.get(QMProperty::wfoverlap);
      U = model->overlaps();
      break;
    }

    case QMProperty::nacvector:{
      arma::mat & nac = props.get(QMProperty::nacvector);
      nac.zeros();

      const arma::uvec &idx = *props.get_idx(QMProperty::nacvector);
      size_t A = idx[0]; size_t B = idx[1];      

      nac(0) = model->nacs()(A,B);
      break;
    }
      
    case QMProperty::qmgradient:{
      arma::mat & g_qm = props.get(QMProperty::qmgradient);
      g_qm.zeros();
      
      arma::uword surface = 0; //assume ground state
      if (props.has_idx(QMProperty::qmgradient)){
	surface = (*props.get_idx(QMProperty::qmgradient))[0];
      }
      g_qm(0) = model->gradients()(surface);

      break;
    }
      
    case QMProperty::qmgradient_multi:{
      arma::cube & g_qm = props.get(QMProperty::qmgradient_multi);
      g_qm.zeros();
      
      arma::uvec surfaces = arma::regspace<arma::uvec>(0,excited_states); //assume all
      if (props.has_idx(QMProperty::qmgradient_multi)){ // specific states
	surfaces = *props.get_idx(QMProperty::qmgradient_multi);
      }

      if (g_qm.n_slices != surfaces.n_elem){
	throw std::range_error("Insufficient space for requested gradients!");
      }
     
      arma::uword i=0;
      for (arma::uword surface: surfaces){
        g_qm.slice(i)(0) = model->gradients()(surface);
        i++;
      }
      break;
    }
      
    case QMProperty::energies:{
      arma::vec & energies = props.get(QMProperty::energies);
      energies.zeros();
      
      if (props.has_idx(QMProperty::energies)){
      	energies = model->energies();
      }
      else{
      	energies = model->energies()(0);
      }
      break;
    }

    case QMProperty::mmgradient: //explicit fall-through
    case QMProperty::mmgradient_multi:
      break;
      
    default:
      throw std::invalid_argument("Unknown QMProperty!");
      break;
    }
  }
}

AvoidedCrossing::AvoidedCrossing(void){
  energy.set_size(2);
  gradient.set_size(2);
  nac.set_size(2,2);
  overlap.set_size(2,2);
}

void  AvoidedCrossing::update(double x){
  x = x+7;
  
  double e = A*std::tanh(B*x);
  double v = C * std::exp(-1.0*x*x);

  energy(1) = std::sqrt(e*e + v*v);
  energy(0) = -1.0 * energy(1);

  double f;
  {
    double sech = 1.0 / std::cosh(B*x);
    f = e*A*B*sech*sech + -2.0*x*v*v;
  }
  gradient = f / energy;
}
