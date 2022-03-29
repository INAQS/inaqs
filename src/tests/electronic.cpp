#include <armadillo>
#include <iostream>
#include "electronic.hpp"
#include <catch.hpp>
#include <complex>
#include "util.hpp"

double TrLg2(arma::mat A){
  return std::real(arma::trace(arma::logmat(A).t()*arma::logmat(A)));
}

TEST_CASE( "Jacobi sweeps locate minimum for A", "[Electronic]"){
  arma::mat A8 = {{-0.6575,  0.3565, -0.6354, -0.1920},
                  {-0.1351, -0.6081, -0.4038,  0.6700},
                  {-0.0916, -0.6991, -0.0847, -0.7041},
                  {-0.7355, -0.1199,  0.6527,  0.1363}};

  arma::vec p = arma::ones(4);
  Electronic::phase_match(A8, p);
  CHECK(TrLg2(A8) == Approx(6.8250));
  REQUIRE((A8*A8.t() - arma::eye(arma::size(A8))).is_zero(1e-14));
}

TEST_CASE( "Jacobi sweeps locate minimum for B", "[Electronic]"){
  arma::mat B8 = {{-0.5987, -0.1138,  0.5219,  0.5969},
                  { 0.5288, -0.5520,  0.6321, -0.1274},
                  { 0.5694,  0.1139, -0.2188,  0.7842},
                  {-0.1942, -0.8182, -0.5294,  0.1121}};

  arma::vec p = arma::ones(4);
  Electronic::phase_match(B8, p);
  CHECK(TrLg2(B8) == Approx(7.7673));
  REQUIRE((B8*B8.t() - arma::eye(arma::size(B8))).is_zero(1e-14));
}

TEST_CASE("More Jacobi Sweeps from ZZ", "[Electronic]"){
  arma::mat T1 = {{-1,  0,  0},
                  { 0,  0, -1},
                  { 0, -1,  0}};
  T1.print("T1 before");

  arma::vec p = arma::ones(3);
  Electronic::phase_match(T1, p);
  T1.print("T1 after");
  CHECK(TrLg2(T1) == Approx(4.9348));

  arma::mat T2 = {{ 0, -1,  0},
                  {-1,  0,  0},
                  { 0,  0, -1}};

  p = arma::ones(3);
  Electronic::phase_match(T2, p);
  CHECK(TrLg2(T2) == Approx(4.9348));
}
