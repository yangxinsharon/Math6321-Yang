// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 4 -- Oct 25
// Problem 2 -- Runge-Kutta-Fehlberg Method
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement: Use both the adaptive forward Euler and adaptive RKF solvers
// to the IVP y'(t) = -y(t) + 2*cos(t), y(0) = 1, t = [0,10]
// for the relative tolerance rtol = {1e-4,1e-6,1e-8} and atol = rtol/100.


#include <iostream>
#include "adapt_euler.hpp"
#include "adapt_rkf.hpp"

using namespace std;
using namespace arma;

class RHS: public RHSFunction {
public:
  int Evaluate(double t, vec& y, vec& f) {    // evaluates the RHS function, f(t,y)
    f = -y + 2.0*cos(t);
    return 0;
  }
};

class RHSt: public RHSFunction {
public:
  int Evaluate(double t, vec& y, vec& ft) {    // evaluates the RHS function, f_t(t,y)
    ft = -2.0*sin(t);
    return 0;
  }
};

class RHSy: public RHSJacobian {
public:
  int Evaluate(double t, vec& y, mat& fy) {    // evaluates the RHS Jacobian, f_y(t,y)
    fy(0) = -1.0;
    return 0;
  }
};

//    Convenience function for analytical solution
vec ytrue(const double t) {
  vec yt(1);
  yt(0) = sin(t) + cos(t);
  return yt;
};


// main routine
int main() {

  // set problem information
  vec y0("1.0");
  double t0 = 0.0;
  double Tf = 10.0;

  // set desired output times
  int Nout = 11;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);

  // create ODE RHS function objects
  RHS  f;
  RHSt ft;
  RHSy fy;

  // create true solution object
  mat Ytrue(1,Nout);
  for (size_t i=0; i<Nout; i++)
    Ytrue.col(i) = ytrue(tspan(i));

  // tolerances
  vec rtols("1.e-4, 1.e-6, 1.e-8");

  // vector for errors at each h
  vec e(rtols.n_elem);


  //---- adaptive forward Euler ----
  cout << "\nAdaptive Forward Euler:\n";

  // create adaptive forward Euler stepper object (will reset rtol and atol before each solve)
  for (size_t ir=0; ir<rtols.size(); ir++) {
    vec atol(1); atol.fill(rtols(ir)/100.0);
    AdaptEuler AE(f, 0.0, atol, y0);
    // set the relative tolerance and absolute tolerance, and call the solver
    cout << "  rtol = " << rtols(ir) << endl;
    cout << "  atol = " << atol(0) << endl;
    AE.rtol = rtols(ir);
    mat Y1 = AE.Evolve(tspan, y0);
    
    // output solution, relative errors, and overall error
    mat Yerr = abs(Y1-Ytrue)/abs(Ytrue);
    e(ir) = Yerr.max();
    cout << "  overall: "
         << "\t steps = " << AE.steps
         << "\t fails = " << AE.fails
         << "\t max relerr = " << e(ir)
         << endl;
  }

  //---- adaptive Runge-Kutta-Fehlberg ----
  cout << "\nAdaptive Runge-Kutta-Fehlberg:\n";
  // create adaptive Runge-Kutta-Fehlberg stepper object (will reset rtol and atol before each solve)
  for (size_t ir=0; ir<rtols.size(); ir++) {
    vec atol(1); atol.fill(rtols(ir)/100.0);
    AdaptRKF AR(f, 0.0, atol, y0);
    // set the relative tolerance and absolute tolerance, and call the solver
    cout << "  rtol = " << rtols(ir) << endl;
    cout << "  atol = " << atol(0) << endl;
    AR.rtol = rtols(ir);
    mat Y2 = AR.Evolve(tspan, y0);
    
    // output solution, errors, and overall error
    mat Yerr = abs(Y2-Ytrue)/abs(Ytrue);
    e(ir) = Yerr.max();
    cout << "  overall: "
         << "\t steps = " << AR.steps
         << "\t fails = " << AR.fails
         << "\t max relerr = " << e(ir)
         << endl;

  }

  // comments on overall results (effciency,accuracy,etc.)
  cout << "\nEffciency:\n";
  cout << "Runge-Kutta-Fehlberg method costs much fewer steps to achieve the desired error than forward euler.\n";
  cout << "It also barely fails comparing to forward euler.\n";

  cout << "\nAccuracy:\n";
  cout << "Based on the same tolerances, Runge-Kutta-Fehlberg achieves slightly smaller errors than forward euler.\n";

  return 0;
}
