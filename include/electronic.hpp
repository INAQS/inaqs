#ifndef __GIFS_ELECTRONIC_HPP
#define __GIFS_ELECTRONIC_HPP

#include <armadillo>

class Electronic{
public:
  Electronic(void) {}
  Electronic(const arma::cx_mat &cs) : amplitudes (cs) {reserve();}
  Electronic(const arma::cx_vec &c)  : amplitudes (c)  {reserve();}
  Electronic(arma::uword NStates, arma::uword NSets=1, arma::uword active=0) {
    reset(NStates, NSets, active);
  }

  void reset(arma::uword NStates, arma::uword NSets, arma::uword active){
    amplitudes.set_size(NStates, NSets);
    amplitudes.zeros();
    amplitudes.row(active).ones();
    reserve();
  }
  
  void set(const arma::cx_mat &cs) {amplitudes = cs;}
  void set(const arma::cx_vec &c) {amplitudes.col(0) = c;}
  
  const std::complex<double> operator()(arma::uword i, arma::uword j) const {return amplitudes(i,j);}
  const std::complex<double> operator()(arma::uword i) const {return amplitudes(i,0);}
  const arma::cx_vec operator()(void) const {return amplitudes.col(0);}
  std::complex<double> & operator[](arma::uword i) {return amplitudes(i,0);}

  const arma::cx_mat get(void) const { return amplitudes; }
  
  void advance(const arma::cx_mat &H, double dt) {advance_exact(H, dt);}
  void advance_rk4(const arma::cx_mat &H, double dt);
  void advance_exact(const arma::cx_mat &H, double dt);

  static void phase_match(arma::mat &U);
  static void phase_match(arma::cx_mat &U);

  static void unitarize(arma::mat &U);
  
private:
  void reserve(void);

  arma::cx_mat amplitudes;
  arma::cx_mat k1, k2, k3, k4;  
};

#endif
