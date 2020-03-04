#include <armadillo>
#include <iostream>
#include "../electronic.hpp"

int main(void){
  // arma::mat A = {{-0.6575,  0.3565, -0.6354, -0.1920},
  // 		 {-0.1351, -0.6081, -0.4038,  0.6700},
  // 		 {-0.0916, -0.6991, -0.0847, -0.7041},
  // 		 {-0.7355, -0.1199,  0.6527,  0.1363}};

  arma_rng::set_seed_random();
  arma::mat A(4,4, arma::fill::randn);
  
  A.print("A");
  //std::cout << "det(A) = " << arma::det(A) << std::endl;
  Electronic::phase_match(A);
  A.print("A'");
  //std::cout << "det(A) = " << arma::det(A) << std::endl;
 }
