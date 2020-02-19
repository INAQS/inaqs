#ifndef __GIFS_ELECTRONIC_HPP
#define __GIFS_ELECTRONIC_HPP

#include <armadillo>

namespace electronic{
  arma::cx_vec propagate_c(arma::cx_vec c, arma::cx_vec diag);
  arma::cx_vec propagate_c(arma::cx_vec c, arma::cx_vec diag, arma::cx_mat off);
  arma::cx_vec rk4_step(arma::cx_mat H, arma::cx_vec c, double dt);
}

#endif
