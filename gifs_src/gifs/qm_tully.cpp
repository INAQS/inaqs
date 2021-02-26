#include <ctime>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <algorithm>
#include "qm_tully.hpp"
#include "properties.hpp"
#include <armadillo>
#include <unordered_map>

ConfigBlockReader QMTully::tully_reader() {
  ConfigBlockReader reader{"tully"};
  reader.add_entry("model", "tully1");
  return reader;
}

QMTully::QMTully(FileHandle& fh, 
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
  ConfigBlockReader reader = tully_reader();
  reader.parse(fh);

  {
    std::string model_in;
    reader.get_data("model", model_in);

    if (model_in == "tully1"){
      model=std::unique_ptr<Tully1>();
    }
    else{
      throw std::runtime_error("Unknown Model: `" + model_in + "`!");
    }
  }
}

void QMTully::update(void){
  QMInterface::update();
  model->update(crd_qm);
}

void QMTully::get_properties(PropMap &props){
  for (QMProperty p: props.keys()){
    switch(p){
      
    case QMProperty::wfoverlap:{
      arma::mat &U = props.get(QMProperty::wfoverlap);
      U = model->overlaps();
      break;
    }

    case QMProperty::nacvector:{
      const arma::uvec &idx = *props.get_idx(QMProperty::nacvector);
      size_t A = idx[0]; size_t B = idx[1];
      arma::mat & nac = props.get(QMProperty::nacvector);
      nac = model->nacs()(A,B);
      break;
    }
      
    case QMProperty::qmgradient:{
      arma::mat & g_qm = props.get(QMProperty::qmgradient);
      arma::uword surface = 0; //assume ground state
      if (props.has_idx(QMProperty::qmgradient)){
	surface = (*props.get_idx(QMProperty::qmgradient))[0];
      }
      g_qm = model->gradients().slice(surface);

      break;
    }
      
    case QMProperty::qmgradient_multi:{
      arma::cube & g_qm = props.get(QMProperty::qmgradient_multi);
      arma::uvec surfaces = arma::regspace<arma::uvec>(0,excited_states); //assume all
      if (props.has_idx(QMProperty::qmgradient_multi)){ // specific states
	surfaces = *props.get_idx(QMProperty::qmgradient_multi);
      }

      if (g_qm.n_slices != surfaces.n_elem){
	throw std::range_error("Insufficient space for requested gradients!");
      }
     
      arma::uword i=0;
      for (arma::uword surface: surfaces){
        g_qm.slice(i) = model->gradients().slice(surface);
        i++;
      }
      break;
    }
      
    case QMProperty::energies:{
      arma::vec & energies = props.get(QMProperty::energies); 
      if (props.has_idx(QMProperty::energies)){
      	energies = model->energies();
      }
      else{
      	energies = model->energies()(0);
      }
      break;
    }

    case QMProperty::mmgradient:
      throw std::invalid_argument("Model Systems do not support QM/MM!");
      break;

    case QMProperty::mmgradient_multi:
      throw std::invalid_argument("Model Systems do not support QM/MM!");
      break;
      
    default:
      throw std::invalid_argument("Unknown QMProperty!");
      break;
    }
  }
}

