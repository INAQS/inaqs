#include "util.hpp"
#include <armadillo>


namespace Util{

  arma::vec center_of_mass(const arma::mat &R, const arma::vec &m){
    // R.B.: sum(A,1) gives the sum of the elements of each row of A
    return arma::sum(R.each_row() % m.t(), 1) / arma::sum(m);
  }

  /*
    For the (3xN) configuration space vectors A, B, compute the sum of
    the cross products of the vectors in each column: Sum(Ai X Bi)
  */
  arma::vec sum_cross(const arma::mat A, const arma::mat B){
    if ((A.n_cols != B.n_cols) || (3 != A.n_rows) || (3 != B.n_rows)){
      throw std::logic_error("Invalid dimension!");
    }

    arma::vec sum(3, arma::fill::zeros);

    for (arma::uword i = 0; i < A.n_cols; i++){
      sum += arma::cross(A.col(i), B.col(i));
    }

    return sum;  
  }
}
