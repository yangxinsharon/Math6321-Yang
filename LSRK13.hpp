/* Low-storage Runge-Kutta（13，4） time stepper class header file.

   Sharon Yang
   Math 6321 @ SMU
   Fall 2020  */

#ifndef LSRK13_DEFINED__
#define LSRK13_DEFINED__

// Inclusions
#include <cmath>
#include "rhs.hpp"

using namespace arma;

// Low-storage Runge-Kutta（13，4） time stepper class
class LSRK13Stepper {

 private:

  RHSFunction *frhs;           // pointer to ODE RHS function
  vec K1, K2, z;               // reused vectors
  mat A, B, c;
 public:

  // number of steps in last call
  unsigned long int nsteps;

  // constructor (sets RHS function pointer, allocates local data)
  LSRK13Stepper(RHSFunction& frhs_, arma::vec& y) {
    frhs = &frhs_;          // set RHSFunction pointer
    nsteps = 0;
    A.load("A13.txt");      // LSRK13 Butcher coefficients
    B.load("B13.txt");
    c.load("c13.txt");

  };

  // Evolve routine (evolves the solution)
  arma::mat Evolve(arma::vec tspan, double h, arma::vec y);

};

#endif
