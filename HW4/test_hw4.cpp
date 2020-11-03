/* Homework 4 testing program.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2020  */

#include <iostream>
#include "adapt_rkf.hpp"

using namespace arma;
using namespace std;


// ODE RHS function class -- instantiates a RHSFunction
//   includes extra routine to evaluate the analytical solution
class ODESystem: public RHSFunction {
public:
  double lambda;
  // evaluates the RHS function, f(t,y)
  int Evaluate(double t, vec& y, vec& f) {
    f(0) = (6.0+4.0*lambda)*y(0) - (4.0+2.0*lambda)*y(1) - (3.0+2.0*lambda)*y(2) + 5.0*exp(-t);
    f(1) = (6.0*lambda-3.0)*y(0) + (2.0-3.0*lambda)*y(1) + (2.0-3.0*lambda)*y(2) + 6.0*exp(-t);
    f(2) = 15.0*y(0) - 10.0*y(1) - 8.0*y(2) + 2.0*exp(-t);
    return 0;
  }
  // computes the true solution, y(t)
  vec TrueSolution(double t) {
    vec yt(3);
    yt(0) = -(9.0+lambda)/(2.0*lambda+2.0)*exp(-t)
          - 4.0*lambda/(lambda+1.0)*exp(lambda*t)
          + 2.5*cos(t) + 0.5*sin(t);
    yt(1) = (lambda-11.0)/(2.0*lambda+2.0)*exp(-t)
          - 6.0*lambda/(lambda+1.0)*exp(lambda*t)
          - 0.5*cos(t) + 2.5*sin(t);
    yt(2) = -1.5*exp(-t) + 5.5*cos(t) - 1.5*sin(t);
    return yt;
  }
};



// main routine
int main() {

  // set up problem
  vec y0("-2.0, -6.0, 4.0");
  double t0 = 0.0;
  double Tf = 5.0;
  double lambda = -20.0;

  // set desired output times
  int Nout = 11;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);

  // create ODE RHS function object
  ODESystem f;
  f.lambda = lambda;

  // create true solution object
  mat Ytrue(3,Nout);
  for (size_t i=0; i<Nout; i++)
    Ytrue.col(i) = f.TrueSolution(tspan(i));

  // tolerances
  vec rtols("1.e-6, 1.e-8, 1.e-10");
  vec atols("1.e-8, 1.e-10, 1.e-12");
  vec atol(3);

  cout << "\nAdaptive RKF solver, steps and errors vs tolerances:\n";
  AdaptRKF RKF(f, 0.0, atol, y0);
  for (int ir=0; ir<rtols.n_elem; ir++) {

    // set the tolerances, and call the solver
    cout << "  rtol = " << rtols(ir) << ",  atol = " << atols(ir);
    RKF.rtol = rtols(ir);
    atol.fill(atols(ir));
    mat Y = RKF.Evolve(tspan, y0);

    // output overall error
    mat Yerr = abs(Y-Ytrue);
    cout << ",  error = " << norm(Yerr,"inf") << ",  steps = "
         << RKF.steps << ",  fails = " << RKF.fails << endl;
  }

  return 0;
}
