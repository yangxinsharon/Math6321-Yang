// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 4 -- Oct 25
// Problem 3 -- Runge-Kutta-Fehlberg and classical Runge-Kutta Methods
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement: Using AdaptRKF solver from the previous problem,
// along with the fourth-order explicit Runge-Kutta time stepper developed 
// in class, do problem 4.12 from the textbook.

#include <iostream>
#include "adapt_rkf.hpp"
#include "adapt_euler.hpp"
#include "erk4.hpp"

using namespace std;
using namespace arma;

class RHS: public RHSFunction {
public:
  double mu;
  double mu_hat;
  int Evaluate(double t, vec& y, vec& f) {
    mu_hat = 1.0 - mu;
    double D1 = pow(pow(y(0)+mu,2) + y(2)*y(2) , 3.0/2.0);
    double D2 = pow(pow(y(0)-mu_hat,2) + y(2)*y(2) , 3.0/2.0);    // evaluates the RHS function, f(t,y)
    f(0) = y(1);
    f(1) = y(0) + 2*y(3) - mu_hat*(y(0)+mu)/D1 - mu*(y(0)-mu_hat)/D2;
    f(2) = y(3);
    f(3) = y(2) - 2*y(1) - mu_hat*(y(2)/D1) - mu*y(2)/D2;
    return 0;
  }
};


// main routine
int main() {

  // time steps to try
  vec h("0.171, 0.0171, 0.00171, 0.000855");

  // set problem information
  vec y0("0.994, 0.0, 0.0, -2.00158510637908252240537862224");
  double t0 = 0.0;
  double Tf = 17.1;

  // set desired output times
  int Nout = 101;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);
  // create ODE RHS function objects
  RHS f;
  f.mu = 0.012277471;

  // tolerances
  double rtol = 1.e-6;
  vec atol("1.e-12");

  //---- adaptive Runge-Kutta-Fehlberg ----
  cout << "\nAdaptive Runge-Kutta-Fehlberg:\n";
  // create adaptive Runge-Kutta-Fehlberg stepper object
  AdaptRKF AR(f, rtol, atol, y0);
  // set the relative tolerance and absolute tolerance, and call the solver
  cout << "  rtol = " << rtol << endl;
  cout << "  atol = " << atol(0) << endl;
  mat Y0 = AR.Evolve(tspan, y0);
  Y0.save("AdaptRKF.txt", raw_ascii);
  cout << "  overall: "
       << "\t steps = " << AR.steps
       << "\t fails = " << AR.fails
       << endl;


  //---- ERK4 ----
  cout << "\nERK4:\n";
  ERK4Stepper ERK4(f, y0);
  cout << "  h = " << h(0) << endl;
  mat Y1 = ERK4.Evolve(tspan, h(0), y0);
  Y1.save("ERK4_100.txt", raw_ascii);

  cout << "  h = " << h(1) << endl;
  mat Y2 = ERK4.Evolve(tspan, h(1), y0);
  Y2.save("ERK4_1000.txt", raw_ascii);

  cout << "  h = " << h(2) << endl;
  mat Y3 = ERK4.Evolve(tspan, h(2), y0);
  Y3.save("ERK4_10000.txt", raw_ascii);

  cout << "  h = " << h(3) << endl;
  mat Y4 = ERK4.Evolve(tspan, h(3), y0);
  Y4.save("ERK4_20000.txt", raw_ascii);

  return 0;
}
