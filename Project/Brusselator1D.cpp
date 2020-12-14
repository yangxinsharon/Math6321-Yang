/* Stability-Optimized Explicit Runge-Kutta Methods:
   3 fourth-order methods: LSRK(12,4), LSRK(13,4), LSRK(14,4) with different stages
   3 test problems: 3) 1D Brusselator
                    
   Sharon Yang
   Math 6321 @ SMU
   Fall 2020  */

#include <iostream>
#include <iomanip>
#include "LSRK12.hpp"
#include "LSRK13.hpp"
#include "LSRK14.hpp"
#include "erk4.hpp"
#include "fwd_euler.hpp"

using namespace std;
using namespace arma;

// ODE RHS function classes
class RHS: public RHSFunction {
public:
  double a = 1.0;
  double b = 3.0;
  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 1.0;
  double k4 = 1.0;  
  int Evaluate(double t, vec& y, vec& f) {    // evaluates the RHS function, f(t,y)
    f(0) = k1*a - k2*b*y(0) + k3*y(0)*y(0)*y(1) - k4*y(0);
    f(1) = k2*b*y(0) - k3*y(0)*y(0)*y(1);
    return 0;
  }
};


// main routine
int main() {

  // time steps to try
  // vec h("0.1 0.05 0.01");
  vec h("0.1");

  // initial condition and time span
  vec y0("1.0 1.0");
  double t0 = 0.0;
  double Tf = 30.0;

  // set desired output times
  int Nout = 301;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);

  // create ODE RHS function object
  RHS f;

  // problem 3
  cout << "\nProblem 3 1D Brusselator:\n";
  char filename[64] = {0};
  //------------LSRK12---------------  
  cout << "\nLSRK12:\n";
  LSRK12Stepper LSRK12(f,y0);
  // create LSRK12 solvers
  for (int ih=0; ih<h.n_elem; ih++) {

    // call stepper
    mat Y = LSRK12.Evolve(tspan, h(ih), y0);
    sprintf(filename, "Brusselator_Y12_%d.txt", ih);
    Y.save(filename,raw_ascii);

  }


  //------------LSRK13---------------
  cout << "\nLSRK13:\n";
  LSRK13Stepper LSRK13(f,y0);
  // create LSRK13 solvers
  for (int ih=0; ih<h.n_elem; ih++) {

    // call stepper
    mat Y = LSRK13.Evolve(tspan, h(ih), y0);

    sprintf(filename, "Brusselator_Y13_%d.txt", ih);
    Y.save(filename,raw_ascii);

  }

  //------------LSRK14---------------
  cout << "\nLSRK14:\n";
  LSRK14Stepper LSRK14(f,y0);
  // create LSRK14 solvers
  for (int ih=0; ih<h.n_elem; ih++) {

    // call stepper
    mat Y = LSRK14.Evolve(tspan, h(ih), y0);

    sprintf(filename, "Brusselator_Y14_%d.txt", ih);
    Y.save(filename,raw_ascii);
  }


  //------------ERK4---------------
  cout << "\nERK4:\n";
  ERK4Stepper ERK4(f,y0);
  // create ERK4 solvers
  for (int ih=0; ih<h.n_elem; ih++) {

    // call stepper

    mat Y = ERK4.Evolve(tspan, h(ih), y0);

    sprintf(filename, "Brusselator_Yerk4_%d.txt", ih);
    Y.save(filename,raw_ascii);
  }


  //------------FE---------------
  cout << "\nFE:\n";
  ForwardEulerStepper FE(f,y0);
  // create FE solvers
  for (int ih=0; ih<h.n_elem; ih++) {

    // call stepper
    mat Y = FE.Evolve(tspan, h(ih), y0);

    sprintf(filename, "Brusselator_Yfe_%d.txt", ih);
    Y.save(filename,raw_ascii);
  }


  tspan.save("brusselator_tspan.txt", raw_ascii);

  return 0;
}
