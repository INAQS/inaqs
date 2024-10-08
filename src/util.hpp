#ifndef __GIFS_UTIL_HPP
#define __GIFS_UTIL_HPP

#include <armadillo>

namespace util{
  arma::vec center_of_mass(const arma::mat &R, const arma::vec &m);
  arma::vec sum_cross(const arma::mat &A, const arma::mat &B);
  //  arma::vec net(arma::vec (*op)(const arma::vec &a, const arma::vec &b), const arma::mat A, const arma::mat B);

  double hypot(std::complex<double> a, std::complex<double> b);
  arma::uword sample_discrete(const arma::vec &p);
  arma::uvec range(arma::uword a, arma::uword b); //[a, b)
  arma::uvec range(arma::uword n); // [0, n)

  arma::mat logmat_unitary(const arma::mat &U);
  arma::cx_mat logmat_unitary(const arma::cx_mat &U);
  void unitarize(arma::mat &U);

  template <typename T>
  bool approx_equal(T a, T b, T tol=arma::datum::eps){
    return std::abs(a-b) < tol;
  };


  /*
    hash functions for comparing arma::mat and friends. Should only be
    compared between containers of the same shape & dimension
  */
  template <template <typename> class C, class T>
  size_t arma_hash(const C<T> & A){
    size_t seed = 0;
    std::hash<T> hash_elem;        // consider element
    for (size_t i = 0; i < A.n_elem; i++){
      seed ^= hash_elem(A(i)) * i; // consider its position
    }
  
    return seed;
  }

}

#endif
