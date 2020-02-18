#include "bomd.hpp"
#include <armadillo>

class FSSH: public BOMD{
public:
  explicit FSSH(int nqm, const int * qmid): BOMD(nqm, qmid) {}; // need to parse our config
  virtual ~FSSH() {};
  
protected:
  void main(void);

  arma::cx_vec c;
  arma::mat U;
  arma::mat T;

  arma::vec V;

  double dtc;
  double dtq;
  
  size_t min_state;
  size_t excited_states;
  size_t active_state;
  size_t target_state;
  bool hopping = false;
  
};
