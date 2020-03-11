#ifndef __GIFS_UTIL_HPP
#define __GIFS_UTIL_HPP

#include <armadillo>

namespace Util{
  arma::vec center_of_mass(const arma::mat R, const arma::vec m);
  arma::vec sum_cross(const arma::mat A, const arma::mat B);
  //  arma::vec net(arma::vec (*op)(const arma::vec &a, const arma::vec &b), const arma::mat A, const arma::mat B);
}

#endif
