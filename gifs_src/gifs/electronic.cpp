#include <armadillo>
#include "electronic.hpp"

namespace electronic{
  // FIXME: maybe should operate directly on c via reference rather than returning? 
  arma::cx_vec rk4_advance(const arma::cx_mat &H, const arma::cx_vec &c, double dt){
    // From NR 17.1.3 and A&S 25.5.10. Require $i \hbar \dot{c} = H c$

    // H_ = -i/hbar * H  for hbar==1 
    arma::cx_mat H_ = H * -1.0 * (std::complex<double> (0, 1));
    arma::cx_vec k1 = dt * H_ * c;
    arma::cx_vec k2 = dt * H_ * (c + (k1 / 2));
    arma::cx_vec k3 = dt * H_ * (c - (k1 / 3) + k2);
    arma::cx_vec k4 = dt * H_ * (c + (k1    ) - k2 + k3);

    return c + (k1 / 6) + (k2 / 3) + (k3 / 3) + (k4 / 6);
  }
}
