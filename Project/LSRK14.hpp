/* Low-storage Runge-Kutta（14，4） time stepper class header file.

   Sharon Yang
   Math 6321 @ SMU
   Fall 2020  */

#ifndef LSRK14_DEFINED__
#define LSRK14_DEFINED__

// Inclusions
#include <cmath>
#include "rhs.hpp"

using namespace arma;

// Low-storage Runge-Kutta（14，4） time stepper class
class LSRK14Stepper {

 private:

  RHSFunction *frhs;           // pointer to ODE RHS function
  vec K1, K2, z;               // reused vectors
  mat A, B, c;
 public:

  // number of steps in last call
  unsigned long int nsteps;

  // constructor (sets RHS function pointer, allocates local data)
  LSRK14Stepper(RHSFunction& frhs_, arma::vec& y) {
    frhs = &frhs_;          // set RHSFunction pointer
    nsteps = 0;
    A.load("A14.txt");  // LSRK14 Butcher coefficients
    B.load("B14.txt");
    c.load("c14.txt");

  };

  // Evolve routine (evolves the solution)
  arma::mat Evolve(arma::vec tspan, double h, arma::vec y);

};

#endif
