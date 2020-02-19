#include <armadillo>
#include "electronic.hpp"

namespace electronic{
  arma::cx_vec rk4_step(arma::cx_mat H, arma::cx_vec c, double dt){
    /*
      From NR 17.1.3 and A&S 25.5.10. Require $\dot{c} = H c$
    */

    arma::cx_vec k1 = dt * H * c;
    arma::cx_vec k2 = dt * H * (c + (k1 / 2));
    arma::cx_vec k3 = dt * H * (c - (k1 / 3) + k2);
    arma::cx_vec k4 = dt * H * (c + (k1    ) - k2 + k3);

    return c + (k1 / 6) + (k2 / 3) + (k3 / 3) + (k4 / 6);
  }
}
