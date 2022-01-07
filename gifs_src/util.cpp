#include "util.hpp"
#include <armadillo>


namespace util{
  /*
    For a discrete probability distribution Pi = p(i), returns the index
    of an element randomly selected by sampling the distribution.

    Examples:
    * sample_discrete({0.5, 0.5}) returns 0 or 1 with probability 0.5
    * sample_discrete({0.0, 0.0, 1.0}) returns 2 with probability 1.0

    Note:
    While there does exist std::discrete_distribution, I think that the
    below implementation is superior in that it does not require
    instantiating a new distribution for every round.
  */
  arma::uword sample_discrete(const arma::vec &p){
    /*
      Sanity checks; require:
      * Pi >= 0 forall i
      * Sum(Pi) == 1
      */
    std::string die = "";
    const double eps = 1e-12;
    if (p.has_nan()){die = "P is malformed; it contains NaN!";}
    if (arma::any(p < -1 * eps)){die = "P cannot have negative elements!";}
    if (std::abs(arma::sum(p) - 1) > eps){ die = "P is not normed!";}

    if (die != ""){
      p.t().print("P");
      throw std::logic_error(die);
    }

    const double zeta = arma::randu();
    return as_scalar(arma::find(arma::cumsum(p) > zeta, 1));
  }


  // return a uvec with filled with the indices [a,b) like Python's
  // range(); contrast arma::span(a,b) -> [a,b]
  arma::uvec range(arma::uword a, arma::uword b){
    if (a > b){
      throw std::logic_error("cannot compute range with a > b");
    }

    arma::uword N = b - a;
    arma::uvec R(N, arma::fill::zeros);
    for (arma::uword i = 0; i < N; i++){
      R(i) = i + a;
    }

    return R;
  }

  // variant for [0, n)
  arma::uvec range(arma::uword n){ return range(0, n); }


  /*
    std::hypot() for complex<double>

    returns sqrt(|a|^2 + |b|^2)

    I'm not sure what is the "best" sequence for the real() operations.
    The present choice---to call real() twice, once after each
    multiplication---was made because that's the first time the
    operands mathematically lie on the real line. Waiting longer seemed
    to invite the accumulation of small numerical deviations.
  */
  double hypot(std::complex<double> a, std::complex<double> b){
    return std::sqrt(std::real(std::conj(a)*a) + std::real(std::conj(b)*b) );
  }


  arma::vec center_of_mass(const arma::mat &R, const arma::vec &m){
    // R.B.: sum(A,1) gives the sum of the elements of each row of A
    return arma::sum(R.each_row() % m.t(), 1) / arma::sum(m);
  }

  /*
    For the (3xN) configuration space vectors A, B, compute the sum of
    the cross products of the vectors in each column: Sum(Ai X Bi)
  */
  arma::vec sum_cross(const arma::mat &A, const arma::mat &B){
    if ((A.n_cols != B.n_cols) || (3 != A.n_rows) || (3 != B.n_rows)){
      throw std::logic_error("Invalid dimension!");
    }

    arma::vec sum(3, arma::fill::zeros);

    for (arma::uword i = 0; i < A.n_cols; i++){
      sum += arma::cross(A.col(i), B.col(i));
    }

    return sum;  
  }


  /*
    Implements Algorithm 3 from:
    Loring, T. A. Computing a logarithm of a unitary matrix with
    general spectrum. Numer. Linear Algebr. Appl. 2014, 21, 744−760

    It is required that the input matrix be approximately unitary.
  */
  arma::mat logmat_unitary(const arma::mat &U){
    return arma::real(logmat_unitary(
      arma::cx_mat(U, arma::zeros(arma::size(U)))
                                     ));
  }


  arma::cx_mat logmat_unitary(const arma::cx_mat &U){
    arma::cx_mat Q,T;
    if ((arma::norm(U.t() * U - arma::eye(arma::size(U))) > 0.75) ||
        !arma::schur(Q, T, U)){
      std::cerr << "Matrix deviates too much from unitary "
        "or Schur Decomposition failed; "
        "falling back to arma::logmat()" << std::endl;
      return arma::logmat(U);
    }

    arma::cx_mat D(arma::size(U), arma::fill::zeros);
    D.diag() = arma::log(T.diag() / arma::abs(T.diag()));
    arma::cx_mat L = Q * D * Q.t();

    return L;
  }


  /*
    Replace U with the closest unitary matrix according to the
    transformation U' = (UU^T)^(-1/2)U.

    For a thorough discussion, see:
    https://en.wikipedia.org/wiki/Kabsch_algorithm
    https://en.wikipedia.org/wiki/Polar_decomposition
  */
  void unitarize(arma::mat &U){
    // arma::cx_mat R; arma::cx_vec s;
    // arma::eig_gen(s, R, U*U.t());
    // U = arma::real(R * diagmat(arma::pow(s, -0.5)) * R.t()) * U;

    /*
      These algorithms are equivalent but the below has better numerical
      stability properties; namely the above returns 2*I for the
      identity (I) rather than I itself.
    */

    arma::mat W,V; arma::vec s;
    arma::svd(W,s,V,U);
    U = W*V.t();
  }
}
