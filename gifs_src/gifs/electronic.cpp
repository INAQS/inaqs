#include <armadillo>
#include "electronic.hpp"


/*
  From NR 17.1.3 and A&S 25.5.10. For systems of the form $i \hbar
  \dot{c} = H c$, advances c using a 4th order Runge-Kutta integrator
*/
void Electronic::advance(const arma::cx_mat &H, double dt){
  k1 = dt * (std::complex<double> (-1, 1)) * H * amplitudes;
  k2 = dt * (std::complex<double> (-1, 1)) * H * (amplitudes + (k1 / 2));
  k3 = dt * (std::complex<double> (-1, 1)) * H * (amplitudes - (k1 / 3) + k2);
  k4 = dt * (std::complex<double> (-1, 1)) * H * (amplitudes + (k1    ) - k2 + k3);

  amplitudes = amplitudes + (k1 / 6) + (k2 / 3) + (k3 / 3) + (k4 / 6);
}

void Electronic::reserve(void){
  auto s = arma::size(amplitudes);
  k1.set_size(s);
  k2.set_size(s);
  k3.set_size(s);
  k4.set_size(s);
}

//FIXME: verify that this isn't returning a copy
arma::cx_vec Electronic::get(void) const{
  return amplitudes.col(0);
}
