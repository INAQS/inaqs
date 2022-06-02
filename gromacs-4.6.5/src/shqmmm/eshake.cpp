#include "eshake.hpp"
#include "gmx_fatal.h"

#include <iostream>
#include <utility>
#include <armadillo>

#include "gifs.hpp"
#include "qm_interface.hpp"
#include "properties.hpp"

#define VERBPRINT(X) do{std::cout << "DVCS [" << __func__ << "]: " << #X << ": " << X << std::endl;}while(0)

void get_econstraint(const arma::mat &R, double & sigma, arma::mat & del_sigma);
void get_econstraint(double & sigma, arma::mat & del_sigma);

double compute_k(arma::mat & V, const arma::rowvec & minv);
double compute_g(arma::mat & R, const arma::rowvec & minv, real invdt);
void bracket_zero(std::function<double(double)> f,
                  const int maxiter, int & niter,
                  double &g0, double & g1, double & f0, double & f1);
void bracket_min(std::function<double(double)> f,
                 const int maxiter, int & niter,
                 double &g0, double & g1, double & f0, double & f1);
void linesearch(std::function<double(double)> f, double thresh_e,
                const int maxiter, int & niter,
                double & g0, double & g1, double & f0, double & f1);

real inaqs_eshake(int natoms, real * invmass_ptr, int econq, real * qs, real invdt){
  if (!gifs_interface_is_ready()){
    gmx_warning("[E-SHAKE] INAQS is not ready; skipping constraint step of type %d", econq);
    return 0;
  }
  
  arma::rowvec minv(natoms, arma::fill::zeros);
  for (int i = 0; i < natoms; i++){minv(i) = invmass_ptr[i];}

  arma::mat Q(3, natoms, arma::fill::zeros);
  for (int i = 0; i < natoms; i++){
    for (int u = 0; u < 3; u++){
      Q(u, i) = qs[i*3+u];
    }
  }

  real lagr;
  switch (econq){
  case (econqCoord):{
    lagr = compute_g(Q, minv, invdt);
    break;
  }
  case (econqVeloc):{
    lagr = compute_k(Q, minv);
    break;
  }
  default:  // should never obtain
    gmx_fatal(FARGS, "[E-SHAKE] Internal error, Electronic SHAKE, called something besieds coords and velocities");
    break;
  }

  // write out finall coords:
  for (int i = 0; i < natoms; i++){
    for (int u = 0; u < 3; u++){
      qs[i*3+u] = Q(u, i);
    }
  }

  return lagr;
}

void get_econstraint(double & sigma, arma::mat & del_sigma){
  arma::cube diabatic_gradients;
  arma::mat diabatic_H;

  PropMap props{};
  // FIXME: how do we know which adiabats? set in inaqs_config.ini?
  props.emplace(QMProperty::diabatic_gradients, &diabatic_gradients);
  props.emplace(QMProperty::diabatic_H, &diabatic_H);
  gifs_QMInterface()->get_properties(props);

  // want to be able to do something like:
  // inaqs::saveh5(diabatic_H)
  // except that we need to make sure we save the right one! not just the last one.
  // similarly, might want to be able to replace gifs_QMInterface() with:
  // inaqs::QMInterface ... 
  
  del_sigma = diabatic_gradients.slice(1) - diabatic_gradients.slice(0);
  sigma = diabatic_H(1,1) - diabatic_H(0,0);
  std::cout << "[E-SHAKE] dE = " << sigma << std::endl;
}

void get_econstraint(const arma::mat &R, double & sigma, arma::mat & del_sigma){
  gifs_update_coords(R);
  get_econstraint(sigma, del_sigma);
}


/*
  Analytic formula for RATTLE velocity correction;
  c.f. DVCS Rattle Notes p. 2
*/
double compute_k(arma::mat & V, const arma::rowvec & minv){
  // FIXME:DVCS: need to ensure properties will be evaluated at *constrained* position; perhaps from last run? perhaps not? looks like it is since the dE rets from k() match the last g()
  // safe: update_gradient(current_position);

  std::cout << "[E-SHAKE] compute_k()" << std::endl;
  
  double sigma;
  arma::mat del_sigma;
  get_econstraint(sigma, del_sigma);

  double vds = arma::as_scalar(V.as_col().t() * del_sigma.as_col());
  double sms = arma::as_scalar(del_sigma.as_col().t() * (minv % del_sigma.each_row()).as_col());

  //FIXME:DVCS tolerance: need to do some thinking about how small vds/sms can be 
  double k = 2 * vds/sms;
  
  V += -k/2*(del_sigma.each_row() % minv);

  double residual = arma::as_scalar(V.as_col().t() * del_sigma.as_col());
  
  std::cout << "[E-SHAKE] k = " << k << "; residual = " << residual << std::endl;
  return  k;
}

/*
  Itterative solution for generalized SHAKE correction;
  c.f. DVCS Rattle Notes p. 2
*/
double compute_g(arma::mat & R, const arma::rowvec & minv, real invdt){
  std::cout << "[E-SHAKE] compute_g()" << std::endl;

  //FIXME:DVCS: want maxiter, thresh, max_step to be configurable
  const int maxiter = 50;  // max # of function evals
  const double thresh_e = 4e-5;  // Ha, ~1mev
  const double max_step = 0.1;   // nm; perhaps should be weighted by sqrt(N)?
  int niter = 0;

  double g = 0;
  double sigma;
  arma::mat del_sigma;
  get_econstraint(sigma, del_sigma); niter++; // FIXME: should be able to pick this up from compute_k()

  //debug
  static int g_call_idx = 0;
  std::cout.precision(11);
  std::cout.setf(std::ios::fixed);

  if (std::abs(sigma) < thresh_e){  // see if we don't need to do anything at all
    std::cout << "[E-SHAKE] " << niter << ": g = " << g << "; sigma = " << sigma << std::endl;
    //R.t().raw_print("Accecpted constrained R for g_call_idx=" + std::to_string(g_call_idx++));
    return g;
  }

  const arma::mat u = R; // 'unconstrained' step
  const arma::mat v = -0.5/invdt * (del_sigma.each_row() % minv);  // search direction
  
  double g0 = g;
  double f0 = sigma;

  // Newton-Raphson guess for bracketing
  static double alpha = 1.0;  // gradients are often small so need alpha 
  double dsdg = arma::as_scalar(v.as_col().t() * del_sigma.as_col());
  g += -alpha * sigma / dsdg;
  std::cout << "[E-SHAKE] " << niter << ": dsdg = " << dsdg << "; a = " << alpha << std::endl;

  {
    double vnorm = arma::norm(v.as_col());
    std::cout << "[E-SHAKE] " << niter << ": |v|=" << vnorm << "; max initial g is: " << max_step/vnorm << "." << std::endl;
    if (vnorm * g > max_step){ // make sure we don't take wild steps
      std::cout << "[E-SHAKE] " << niter << ": guess step too large, " << g*vnorm << "; resetting to " << max_step << "." << std::endl;
      g = max_step / vnorm;
    }
  }
  
  get_econstraint(u + g*v, sigma, del_sigma); niter++; // FIXME: don't need this gradient eval
  double g1 = g;
  double f1 = sigma;
  
  if (std::abs(sigma) < thresh_e){ // maybe our first guess is already good enough!
    std::cout << "[E-SHAKE] " << niter << ": g = " << g << "; sigma = " << sigma << std::endl;
    R = u + g*v;
    //R.t().raw_print("Accecpted constrained R for g_call_idx=" + std::to_string(g_call_idx++));
    return g;
  }

  // auto-scaling alpha over the course of the run
  const double alpha_fac = 1.2; // will walk [1e0, 1e-1] in 12 steps
  if (f0*f1 > 0 && alpha < 1.0/alpha_fac){  // if we didn't go far enough and have headroom, stretch!
    alpha *= alpha_fac;
  }
  else if (std::abs(f1) > thresh_e){ // shrink as long as we're not at thresh yet
    alpha /= alpha_fac;
  }
  
  // function for bracketing and line search
  auto ediff = [u,v](double g) -> double{
    arma::mat temp;  // FIXME: Don't need gradient so don't compute!
    double s;
    get_econstraint(u+g*v, s, temp);
    return s;
  };
  
  bracket_zero(ediff, maxiter, niter, g0, g1, f0, f1);
  //bracket_min(ediff, maxiter, niter, g0, g1, f0, f1);  // for when we want to do CI seam
  linesearch(ediff, thresh_e, maxiter, niter, g0, g1, f0, f1);

  // Take the smallest one:
  if (std::abs(f0) < std::abs(f1)){
    sigma = f0;
    g = g0;
  }
  else{
    sigma = f1;
    g = g1;
  }
  
  std::cout << "[E-SHAKE] " << niter << ": g = " << g << "; sigma = " << sigma << std::endl;
  R = u + g*v;


  //R.t().raw_print("Accecpted constrained R for g_call_idx=" + std::to_string(g_call_idx++));
  
  return g;
}


// Straight out of Numerical Recipes, Press et al, sec 9.1
void bracket_zero(std::function<double(double)> f,
                  const int maxiter, int & niter,
                  double &g0, double & g1, double & f0, double & f1){
  const double FACTOR = 1.6;

  if (f1<f0){ // make sure we are (low, high)
    std::swap(g1,g0);
    std::swap(f1,f0);
  }

  std::cout << "[E-SHAKE] " << niter << ": ";
  std::cout << "guess: g is bracketed by (" << g0 << "," << g1 << ") ";
  std::cout << "with f=(" << f0 << "," << f1 << ")" << std::endl;
  
  while (f0*f1 > 0){
    if (std::abs(f0) < std::abs(f1)){
      g0 += FACTOR * (g0-g1);
      f0 = f(g0); niter++;
    }
    else{
      g1 += FACTOR * (g1-g0);
      f1 = f(g1); niter++;
    }
    std::cout << "[E-SHAKE] " << niter << ": ";
    std::cout << "bracketing: g is bracketed by (" << g0 << "," << g1 << ") ";
    std::cout << "with f=(" << f0 << "," << f1 << ")" << std::endl;

    if (niter > maxiter){
      gmx_fatal(FARGS, "[E-SHAKE] bracket failes to converge!");
    }
  }
}

// Straight out of Numerical Recipes, Press et al, sec 10.1
// FIXME:DVCS need to write this so we can sample CI seams too
void bracket_min(std::function<double(double)> f,
                 const int maxiter, int & niter,
                 double &g0, double & g1, double & f0, double & f1){

}

// linesearch with false-position (bisection with assumption of linearity)
void linesearch(std::function<double(double)> f, double thresh_e,
                  const int maxiter, int & niter,
                  double & g0, double & g1, double & f0, double & f1){
  double sigma, g;
  
  if (std::abs(f0) < std::abs(f1)){
    sigma = f0;
    g = g0;
  }
  else{
    sigma = f1;
    g = g1;
  }

  double dg = g1 - g0;
  while (std::abs(sigma) > thresh_e){
    //g = g0 + (dg *= 0.5); // bisection
    g = g0 + (dg = g1 - g0)*f0/(f0-f1); // false-position; seems much faster
    sigma = f(g); niter++;
    if (sigma < 0) {g0 = g; f0 = sigma;}
    else           {g1 = g; f1 = sigma;}

    std::cout << "[E-SHAKE] " << niter << ": ";
    std::cout << "bisecting: g is bracketed by (" << g0 << "," << g1 << ") ";
    std::cout << "with f=(" << f0 << "," << f1 << ") ";
    std::cout << "last dg=" << dg << std::endl;

    if (std::abs(dg) < 1e-10){
      std::cout << "[E-SHAKE] " << niter << ": ";
      std::cout << "Terminating bisection for tiny dg." << std::endl;
      break;
    }
    
    if (niter > maxiter){
      gmx_fatal(FARGS, "[E-SHAKE] bisection failes to converge!");
    }
  }
}
