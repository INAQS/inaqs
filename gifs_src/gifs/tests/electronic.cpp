#include <armadillo>
#include <iostream>
#include "electronic.hpp"
#include <catch.hpp>
#include <complex>

TEST_CASE( "2: A8", "[multi-file:2]"){
  arma::mat A8 = {{-0.6575,  0.3565, -0.6354, -0.1920},
                  {-0.1351, -0.6081, -0.4038,  0.6700},
                  {-0.0916, -0.6991, -0.0847, -0.7041},
                  {-0.7355, -0.1199,  0.6527,  0.1363}};
  Electronic::phase_match(A8);
  REQUIRE(std::real(arma::trace(arma::logmat(A8).t()*arma::logmat(A8))) < 8);
}

// int main(void){

//   std::cout << "Tr[|log(A8)|^2]: " << arma::trace(arma::logmat(A8).t()*arma::logmat(A8)) << std::endl;
//   std::cout << "|A8|: " << arma::det(A8) << std::endl;
//   A8.print("A8");
//   //FIXME: should implement phase_match w/out unitarize s.t. we can do matched = phase_match(unitarize(A));
//   electronic::phase_match(A8);
//   std::cout << "Tr[|log(A8)|^2]: " << arma::trace(arma::logmat(A8).t()*arma::logmat(A8)) << std::endl;
//   std::cout << "|A8|: " << arma::det(A8) << std::endl;
//   A8.print("A1?");

//   arma::mat I = A8 * A8.t();
//   I.clean(1e-15);
//   I.print("I?");

//   return 0;
// }
