#include <armadillo>
#include "electronic.hpp"
#include "util.hpp"

/*
  From NR 17.1.3 and A&S 25.5.10. For systems of the form $i \hbar
  \dot{c} = H c$, advances c using a 4th order Runge-Kutta integrator
*/
// FIXME: Returns different results than exact propagation
void Electronic::advance_rk4(const arma::cx_mat &H, double dt){
  const std::complex<double> midt(0,-dt);

  k1 = midt * H * amplitudes;
  k2 = midt * H * (amplitudes + (k1 / 2));
  k3 = midt * H * (amplitudes - (k1 / 3) + k2);
  k4 = midt * H * (amplitudes + (k1    ) - k2 + k3);

  amplitudes = amplitudes + (k1 / 6) + (k2 / 3) + (k3 / 3) + (k4 / 6);
}


void Electronic::advance_exact(const arma::cx_mat &H, double dt){
  const std::complex<double> midt(0,-dt);
  amplitudes = arma::expmat(midt * H) * amplitudes;
}


void Electronic::reserve(void){
  auto s = arma::size(amplitudes);
  k1.set_size(s);
  k2.set_size(s);
  k3.set_size(s);
  k4.set_size(s);
}


/*
  Implements Zeyu Zhou et al. JCTC 2020, 16, 835--846

  Approximately minimizes Tr[log(U)^2] via Jacobi sweeps while
  enforcing det(U) == 1

  Results agree with the matrices in the paper's Appendix D
*/
void Electronic::phase_match(arma::mat &U, arma::vec &phases){
  U = arma::diagmat(phases) * U; // update phases from previous time step if provided
  phases.ones();// assume we're aligned this time

  // align and construct phases for next step
  for (arma::uword i = 0; i < U.n_cols; i++){
    if (U(i,i) < 0){
      U.col(i)*= -1.0;
      phases(i) *= -1.0;
    }
  }

  // Step 1: det(U) == 1
  util::unitarize(U);
  if (arma::det(U) < 0){
    U.col(0) *= -1.0;
    phases(0) *= -1.0;
  }

  // Step 2: Jacobi sweeps
  double deljk, Ujj, Ukk;

  const arma::uword N = U.n_rows;
  bool change;
  do{
    change = false;
    for (arma::uword j = 0; j < N; j++){
      for (arma::uword k = j + 1; k < N; k++){
	Ujj = U(j,j);
	Ukk = U(k,k);
	// Eq. 34 with the sum implemented as a dot product
	deljk =
            3*(Ujj*Ujj + Ukk*Ukk)
	  + 6*(U(j,k)*U(k,j))
	  + 8*(Ujj + Ukk)
	  - 3*(arma::as_scalar(U.row(j)*U.col(j) + U.row(k)*U.col(k)));
	if (deljk < 0){
	  U.col(j) *= -1.0;
	  U.col(k) *= -1.0;
          phases(j) *= -1.0;
          phases(k) *= -1.0;
          change = true;
	}
      }
    }
  }while(change);
}

void Electronic::phase_match(arma::cx_mat &U, arma::cx_vec &phases){
  (void) U;
  (void) phases;
  throw std::runtime_error("complex phase matching not implemented");
}
