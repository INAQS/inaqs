#ifndef GIFS_QM_MODEL_HPP
#define GIFS_QM_MODEL_HPP

#include "properties.hpp"
#include "qm_interface.hpp"
#include "configreader.hpp"
#include <armadillo>
#include <memory>

class HamiltonianDynamics{
public:
  HamiltonianDynamics(void) {}
  virtual ~HamiltonianDynamics() {};

  virtual void update(double x) =0;
  void scan(double a, double b, int N);
  
  const arma::vec & energies(void)  const {return energy;}
  const arma::vec & gradients(void) const {return gradient;}
  const arma::mat & nacs(void)      const {return nac;}
  const arma::mat & overlaps(void)  const {return overlap;}

protected:
  // updates overlap and, phase matches the new eigen vectors to the old ones.
  void update_overlap(arma::mat &evec);
  void save_state(double x);
  
  arma::vec energy;
  arma::vec gradient;
  arma::mat nac;
  arma::mat overlap;

  bool record = true;
  std::string state_file = "./out.dat";
};


/*
  Implments system A, eq. 1, from Subotnik and Shenvi (2011)
*/
class AvoidedCrossing: public HamiltonianDynamics{
public:
  AvoidedCrossing(void);
  void update(double x) override;
    
private:
  const double A = 0.03;
  const double B = 1.6;
  const double C = 0.005;
};


/*
  Implments system B, eq. 3, from Subotnik and Shenvi (2011)
*/
class ReflectiveAvoidedCrossing: public HamiltonianDynamics{
public:
  ReflectiveAvoidedCrossing(void);
  void update(double x) override;

private:
  const double A = 0.03;
  const double B = 1.6;
  double C = 0.005;
};



class QM_Model: public QMInterface{
public:
  QM_Model(FileHandle& fh, 
           const arma::uvec& in_qmids, 
	   arma::mat& in_qm_crd, 
	   arma::mat& in_mm_crd, 
	   arma::vec& in_mm_chg, 
	   const int charge, 
	   const int mult,
	   const int excited_states,
           const int min_state);

  ~QM_Model(){ delete model; }

  
  void get_properties(PropMap &props) override;
  void update(void) override;

private:
  ConfigBlockReader model_reader();
  HamiltonianDynamics * model;
};

#endif
