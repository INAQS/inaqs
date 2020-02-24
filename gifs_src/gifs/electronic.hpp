#ifndef __GIFS_ELECTRONIC_HPP
#define __GIFS_ELECTRONIC_HPP

#include <armadillo>

namespace electronic{
  /*
    For systems of the form $i \hbar \dot{c} = H c$, advances c using
    a 4th order Runge-Kutta integrator
  */
  arma::cx_vec rk4_advance(const arma::cx_mat &H, const arma::cx_vec &c, double dt);
}

#endif
