/* Adaptive Runge-Kutta-Fehlberg solver class implementation file.

   Class to perform adaptive time evolution of the IVP
        y' = f(t,y),  t in [t0, Tf],  y(t0) = y0
   using the Runge-Kutta-Fehlberg time stepping method.

   Sharon Yang
   Math 6321 @ SMU
   Fall 2020  */

#include "adapt_rkf.hpp"

using namespace std;
using namespace arma;

// Adaptive Runge-Kutta-Fehlberg class constructor routine
//
// Inputs:  frhs_ holds the ODE RHSFunction object, f(t,y)
//          rtol holds the desired relative solution accuracy
//          atol holds the desired absolute solution accuracy
//
// Sets default values for adaptivity parameters, all of which may
// be modified by the user after the solver object has been created
AdaptRKF::AdaptRKF(RHSFunction& frhs_, double rtol_,
                       arma::vec& atol_, arma::vec& y) {
  frhs = &frhs_;    // set RHSFunction pointer
  atol = &atol_;    // set absolute tolerance pointer
  rtol = rtol_;     // copy relative tolerance
  fn = vec(y);      // clone y to create local vectors
  y1 = vec(y);
  y2 = vec(y);
  yerr = vec(y);
  w = vec(y);

  maxit = 1e6;      // set default solver parameters
  bias = 2.0;
  grow = 50.0;
  safe = 0.95;
  ONEMSM = 1.0 - sqrt(eps(1.0));
  ONEPSM = 1.0 + sqrt(eps(1.0));
  p = 5;
  fails = 0;
  steps = 0;
  error_norm = 0.0;
  h = 0.0;

  k0 = vec(y);  // allocate reusable data
  k1 = vec(y);
  k2 = vec(y);
  k3 = vec(y);
  k4 = vec(y);
  k5 = vec(y);
  z = vec(y);   //   based on size of y
  A = mat(6,6); A.fill(0.0); 
  A(1,0) = 1.0/4.0;                
  A(2,0) = 3.0/32.0;
  A(2,1) = 9.0/32.0;
  A(3,0) = 1932.0/2197.0;                
  A(3,1) = -7200.0/2197.0;
  A(3,2) = 7296.0/2197.0;                
  A(4,0) = 439.0/216.0;          
  A(4,1) = -8.0;
  A(4,2) = 3680.0/513.0;
  A(4,3) = -845.0/4104.0;          
  A(5,0) = -8.0/27.0;
  A(5,1) = 2.0;
  A(5,2) = -3544.0/2565.0;          
  A(5,3) = 1859.0/4104.0;
  A(5,4) = -11.0/40.0;
  c = vec(6);
  c(0) = 0.0;
  c(1) = 1.0/4.0;
  c(2) = 3.0/8.0;
  c(3) = 12.0/13.0;
  c(4) = 1.0;
  c(5) = 1.0/2.0;
  b1 = vec(6);
  b1(0) = 25.0/216.0;
  b1(1) = 0.0;
  b1(2) = 1408.0/2565.0;
  b1(3) = 2197.0/4104.0;
  b1(4) = -1.0/5.0;
  b1(5) = 0.0;
  b2 = vec(6);
  b2(0) = 16.0/135.0;
  b2(1) = 0.0;
  b2(2) = 6656.0/12825.0;
  b2(3) = 28561.0/56430.0;
  b2(4) = -9.0/50.0;
  b2(5) = 2.0/55.0;  
};

// Error weight vector utility routine
void AdaptRKF::error_weight(vec& y, vec& w) {
  for (size_t i=0; i<y.size(); i++)
    w(i) = bias / ((*atol)(i) + rtol * std::abs(y(i)));
}


// The adaptive Runge-Kutta-Fehlberg time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          y holds the initial condition, y(t0)
// Outputs: the output matrix holds the computed solution at
//          all tspan values,
//            [y(t0), y(t1), ..., y(tN)]
mat AdaptRKF::Evolve(vec tspan, vec y) {

  // store sizes
  size_t m = y.n_elem;
  size_t N = tspan.n_elem-1;

  // initialize output
  mat Y(m, N+1, fill::zeros);
  Y.col(0) = y;

  // reset counters, current time value
  fails = steps = 0;
  double t = tspan(0);

  // check for legal time span
  for (size_t tstep=0; tstep<N; tstep++) {
    if (tspan(tstep+1) < tspan(tstep)) {
      cerr << "AdaptRKF::Evolve Illegal tspan\n";
      return Y;
    }
  }

  // initialize error weight vector, and check for legal tolerances
  error_weight(y,w);

  // get ||y'(t0)||
  if (frhs->Evaluate(t, y, fn) != 0) {
    std::cerr << "Evolve error in RHS function\n";
    return Y;
  }

  // estimate initial h value via linearization, safety factor
  error_norm = norm(fn%w,"inf");
  error_norm = std::max(error_norm, 1.e-8);
  h = safe/error_norm;

  // iterate over output times
  for (size_t tstep=0; tstep<N; tstep++) {

    // loop over internal steps to reach desired output time
    while ((tspan(tstep+1)-t) > sqrt(eps(tspan(tstep+1)))) {

      // enforce maxit
      if (steps+fails > maxit) {
        std::cerr << "Exceeded maximum allowed iterations\n";
        return Y;
      }

      // bound internal time step
      h = std::min(h,tspan(tstep+1)-t);

      // reset both solution approximations to current solution
      y1 = y;
      y2 = y;

      // stage 0: compute k0 at old step solution
      if (frhs->Evaluate(t, y, k0) != 0) {
        std::cerr << "AdaptiveRKF::Evolve: Error in ODE RHS function\n";
        return Y;
      }

      // stage 1: compute z1 (store in z) and corresponding k1
      z = y + (A(1,0)*h)*k0;
      if (frhs->Evaluate(t+c(1)*h, z, k1) != 0) {
        std::cerr << "AdaptiveRKF::Evolve: Error in ODE RHS function\n";
        return Y;
      }

      // stage 2: compute z2 (store in z) and corresponding k2
      z = y + (A(2,0)*h)*k0 + (A(2,1)*h)*k1;
      if (frhs->Evaluate(t+c(2)*h, z, k2) != 0) {
        std::cerr << "AdaptiveRKF::Evolve: Error in ODE RHS function\n";
        return Y;
      }

      // stage 3: compute z3 (store in z) and corresponding k3
      z = y + (A(3,0)*h)*k0 + (A(3,1)*h)*k1 + (A(3,2)*h)*k2;
      if (frhs->Evaluate(t+c(3)*h, z, k3) != 0) {
        std::cerr << "AdaptiveRKF::Evolve: Error in ODE RHS function\n";
        return Y;
      }

      // stage 4: compute z4 (store in z) and corresponding k4
      z = y + (A(4,0)*h)*k0 + (A(4,1)*h)*k1 + (A(4,2)*h)*k2 + (A(4,3)*h)*k3;
      if (frhs->Evaluate(t+c(4)*h, z, k4) != 0) {
        std::cerr << "AdaptiveRKF::Evolve: Error in ODE RHS function\n";
        return Y;
      }

      // stage 5: compute z5 (store in z) and corresponding k5
      z = y + (A(5,0)*h)*k0 + (A(5,1)*h)*k1 + (A(5,2)*h)*k2 + (A(5,3)*h)*k3 + (A(5,4)*h)*k4;
      if (frhs->Evaluate(t+c(5)*h, z, k5) != 0) {
        std::cerr << "AdaptiveRKF::Evolve: Error in ODE RHS function\n";
        return Y;
      }

      y1 += h*(b1(0)*k0 + b1(1)*k1 + b1(2)*k2 + b1(3)*k3 + b1(4)*k4 + b1(5)*k5);
      y2 += h*(b2(0)*k0 + b2(1)*k1 + b2(2)*k2 + b2(3)*k3 + b2(4)*k4 + b2(5)*k5);

      // compute error estimate
      yerr = y2 - y1;

      // compute error estimate success factor
      error_norm = norm(yerr%w,"inf");
      error_norm = std::max(error_norm, 1.e-8);

      // check error
      if (error_norm < ONEPSM) {  // successful step

        // update current time, solution, error weights, and work counter
        t += h;
        y = y2;
        error_weight(y,w);
        steps++;

      } else {                    // failed step
        fails++;
      }

      // pick next time step size based on this error estimate
      double eta = safe*std::pow(error_norm, -1.0/(p+1));   // step size growth factor
      eta = std::min(eta, grow);                            // limit maximum growth
      h *= eta;                                             // update h

    }

    // store updated solution in output array
    Y.col(tstep+1) = y;

  }

  return Y;
}
