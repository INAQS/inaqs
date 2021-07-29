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
    else if (model_in == "reflective-avoided-crossing"){
      model= new ReflectiveAvoidedCrossing();
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


// updates overlap and, phase matches the new eigen vectors to the old ones.
void HamiltonianDynamics::update_overlap(arma::mat &evec){
  static arma::mat evec_last = evec;
  overlap = evec_last.t() * evec;

  while (arma::any(overlap.diag() < 0)){ // phase match
    arma::uword c = overlap.diag().index_min();
    evec.col(c) *= -1.0;
    
    overlap = evec_last.t() * evec;
  }
  
  evec_last = evec;
}

void HamiltonianDynamics::save_state(double x){
  std::vector<double> state;
     
  state.push_back( x );
  state.push_back( energy(0) );
  state.push_back( energy(1) );
  state.push_back( gradient(0) );
  state.push_back( gradient(1) );
  state.push_back( overlap(0,1) );
  state.push_back( nac(0,1) );
    
  std::ofstream stream;
  stream.open(state_file, std::ios::out | std::ios::app | std::ios::binary);
  if (!stream){
    throw std::runtime_error("Cannot open data file!");
  }
  else{
    for (auto e: state){
      stream << e << " ";
    }
    stream << std::endl;
  }
}

void HamiltonianDynamics::scan(double a, double b, int N){
  record = true;

  if ((b-a < 0) || (N < 0)){
    throw::std::logic_error("Can only scan on [a,b) in N>0 steps!");
  }
  
  const double dx = (b-a)/N;

  for (double x = a; x < b; x += dx){
    update(x);
  }
  
}
  
ReflectiveAvoidedCrossing::ReflectiveAvoidedCrossing(void){
  energy.set_size(2);
  gradient.set_size(2);
  nac.set_size(2,2);
  overlap.set_size(2,2);
}

void ReflectiveAvoidedCrossing::update(double x){
  double e2 = A*std::tanh(B*(x+7));
  double e1 = e2 + 2*A*std::tanh(B*x) + 2*A;
  
  double v = C * std::exp(-1.0*(x+7)*(x+7));

  arma::mat H = {{e1, v },
                 {v, -e2}};

  arma::vec eval;
  arma::mat evec;
  eig_sym(eval, evec, H);

  update_overlap(evec);

  energy = eval;

  // compute gradients
  {
    double de2 = A*B / std::pow(std::cosh(B*(x+7)),2) ;
    double de1 = A*B / std::pow(std::cosh(B*x)    ,2) + de2 ;
    double dv = -C * 2 * (x+7) *  std::exp(-1.0*(x+7)*(x+7));

    // gradient = 0.5*((de1 - de2) +
    //                 ((de2-de1)*(e2-e1) + 2*(e1*de2 + e2*de1) + 4*v*dv) /
    //                 ((e2-e1) + 2 * energy));

    
    double gl = 0.5 * (de1 - de2);
    double gr = 0.5 * ((de2-de1)*(e2-e1) + 2*(e1*de2 + e2*de1) + 4*v*dv) /
      std::sqrt((e2-e1)*(e2-e1) + 4*e1*e2+4*v*v);

    gradient(0) = gl - gr;
    gradient(1) = gl + gr;
  }


  // compute NAC
  nac = evec.t() * arma::diagmat(gradient) * evec;
  nac.diag().zeros();
  {
    double gap = energy(1) - energy(0);
    nac(0,1) /=  gap;
    nac(1,0) /= -gap;
  }
  

  if(record){
    save_state(x);
  }
}



AvoidedCrossing::AvoidedCrossing(void){
  energy.set_size(2);
  gradient.set_size(2);
  nac.set_size(2,2);
  overlap.set_size(2,2);
}

void AvoidedCrossing::update(double x){
  double xs = x+7;
  
  double e = A * std::tanh(B*xs);
  double v = C * std::exp(-1.0*xs*xs);

  arma::mat H = {{e,  v},
                 {v, -e}};

  arma::vec eval;
  arma::mat evec;
  eig_sym(eval, evec, H);

  update_overlap(evec);

  energy = eval;

  // compute gradients
  double f;
  {
    double sech = 1.0 / std::cosh(B*xs);
    f = e*A*B*sech*sech + -2.0*xs*v*v;
  }
  gradient = f / energy;


  // compute NAC
  nac = evec.t() * arma::diagmat(gradient) * evec;
  nac.diag().zeros();
  {
    double gap = energy(1) - energy(0);
    nac(0,1) /=  gap;
    nac(1,0) /= -gap;
  }

  // save a bunch of things
  if(record)
  {
    save_state(x);
  }
}
