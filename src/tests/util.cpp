#include <armadillo>
#include <iostream>
#include <catch.hpp>
#include "util.hpp"


TEST_CASE( "Random unitary matrix logarithms", "[Util]"){
  for (arma::uword N = 2; N < 100; N*=1.5){
    for (int i = 0; i < 10; i++){
      arma::mat A(N,N, arma::fill::randu);
      util::unitarize(A);

      CHECK(arma::norm(util::logmat_unitary(A) -
                       arma::real(arma::logmat(A)))
            == Approx(0.0).margin(1e-9));
    }
  }
}
