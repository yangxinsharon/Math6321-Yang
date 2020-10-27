/* Adaptive Runge-Kutta-Fehlberg time stepper class header file.

   Sharon Yang
   Math 6321 @ SMU
   Fall 2020  */

#ifndef ADAPT_RKF_DEFINED__
#define ADAPT_RKF_DEFINED__

// Inclusions
#include <cmath>
#include "rhs.hpp"

// Adaptive Runge-Kutta-Fehlberg time stepper class
class AdaptRKF {

 private:

  // private data
  RHSFunction *frhs;       // pointer to ODE RHS function
  arma::vec k0, k1, k2, k3, k4, k5, z; // reused vectors
  arma::mat A;             // Butcher table
  arma::vec b1, b2, c;     
  
  arma::vec *atol;         // pointer to user-provided absolute tolerance array
  arma::vec fn;            // local vector storage
  arma::vec y1;
  arma::vec y2;
  arma::vec yerr;
  arma::vec w;
  double grow;        // maximum step size growth factor
  double bias;        // bias factor for error estimate
  double safe;        // safety factor for step size estimate
  double ONEMSM;      // safety factors for
  double ONEPSM;      //    floating-point comparisons
  int p;              // order of accuracy for the method

 public:

  // publicly-accesible input parameters
  double rtol;        // desired relative solution error
  long int maxit;     // maximum allowed number of steps

  // outputs
  double error_norm;  // current estimate of the local error ratio
  double h;           // current time step size
  long int fails;     // number of failed steps
  long int steps;     // number of successful steps

  // constructor (sets RHS function pointer & solver parameters, copies y for local data)
  AdaptRKF(RHSFunction& frhs_, double rtol_, arma::vec& atol_, arma::vec& y);

  // utility routine to compute the error weight vector
  void error_weight(arma::vec& y, arma::vec& w);

  // Evolve routine (evolves the solution via adaptive Runge-Kutta Fehlberg)
  arma::mat Evolve(arma::vec tspan, arma::vec y);

};

#endif
