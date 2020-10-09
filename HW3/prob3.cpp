// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 3 -- Oct 09
// Problem 3 -- Implicit Methods
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement: Solve the problem 
// y' = lambda*y + 1/(1+pow(t,2)) - lambda*arctan(t), y(0) = 0,
// that has analytical solution y(t) = arctan(t).
// Use stiffness parameters lambda = {-2000,-20000,-200000},
// step sizes h = {0.1,0.05,0.025,0.0125}, and
// time interval t = [0,1] for
// (a) generalized midpoint with theta = 1,
// (b) generalized midpoint with theta = 0.5, and 
// (c) generalized midpoint with theta = 0.45,

#include <iostream>
#include "GeneralizedMidpoint.hpp"

using namespace std;
using namespace arma;


// Define classes to compute the ODE RHS function and its Jacobian

//    ODE RHS function class -- instantiates a RHSFunction
class MyRHS: public RHSFunction {
public:
  double lambda;                              // stores some local data
  int Evaluate(double t, vec& y, vec& f) {    // evaluates the RHS function, f(t,y)
    f(0) = lambda*y(0) + 1.0/(1.0+t*t) - lambda*atan(t);
    return 0;
  }
};

//    ODE RHS Jacobian function class -- instantiates a RHSJacobian
class MyJac: public RHSJacobian {
public:
  double lambda;                              // stores some local data
  int Evaluate(double t, vec& y, mat& J) {    // evaluates the RHS Jacobian, J(t,y)
    J(0,0) = lambda;
    return 0;
  }
};


// Convenience function for analytical solution
vec ytrue(const double t) {
  vec yt(1);
  yt(0) = atan(t);
  return yt;
};



// main routine
int main() {

  // time steps to try
  vec h("0.1, 0.05, 0.025, 0.0125");

  // lambda values to try
  vec lambdas("-2000.0, -20000.0, -200000.0");

  // set problem information
  vec y0("0.0");
  double t0 = 0.0;
  double Tf = 1.0;

  // set desired output times
  int Nout = 11;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);

  // create ODE RHS and Jacobian objects
  MyRHS rhs;
  MyJac Jac;

  // create true solution results
  mat Ytrue(1,Nout);
  for (size_t i=0; i<Nout; i++) {
    Ytrue.col(i) = ytrue(tspan(i));
  }


  //------ run first test: theta = 1.0 ------
  double theta = 1.0;
  // create time stepper object
  GeneralizedMidpoint GM(rhs, Jac, theta, y0);

  // update Newton solver parameters
  GM.newt.tol = 1.e-3;
  GM.newt.maxit = 20;
  GM.newt.show_iterates = false;

  // loop over lambda values
  for (int il=0; il<lambdas.n_elem; il++) {

    // set current lambda value into rhs and Jac objects
    rhs.lambda = lambdas(il);
    Jac.lambda = lambdas(il);

    // loop over time step sizes
    for (int ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "\nRunning generalized midpoint with stepsize h = " << h(ih)
	   << ", lambda = " << lambdas(il) << ":\n";
      mat Y = GM.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      /*for (size_t i=0; i<Nout; i++)
        printf("    y(%.1f) = %9.6f   |error| = %.2e\n",
               tspan(i), Y(0,i), Yerr(0,i));*/
      Y.save("Y1.txt", raw_ascii);
      cout << "  Max error = " << Yerr.max() << endl;
      cout << "  Total internal steps = " << GM.nsteps << endl;
      cout << "  Total Newton iterations = " << GM.nnewt << endl;
    }
    cout << endl;
  }


/*
  //------ run second test: theta = 0.50 ------
  double theta = 0.5;
  // create time stepper object
  GeneralizedMidpoint GM(rhs, Jac, theta, y0);

  // update Newton solver parameters
  GM.newt.tol = 1.e-3;
  GM.newt.maxit = 20;
  GM.newt.show_iterates = false;
  // loop over lambda values
  for (int il=0; il<lambdas.n_elem; il++) {

    // set current lambda value into rhs and Jac objects
    rhs.lambda = lambdas(il);
    Jac.lambda = lambdas(il);

    // loop over time step sizes
    for (int ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "\nRunning generalized midpoint with stepsize h = " << h(ih)
     << ", lambda = " << lambdas(il) << ":\n";
      mat Y = GM.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      for (size_t i=0; i<Nout; i++)
        printf("    y(%.1f) = %9.6f   |error| = %.2e\n",
               tspan(i), Y(0,i), Yerr(0,i));
      Y.save("Y2.txt", raw_ascii);
      cout << "  Max error = " << Yerr.max() << endl;
      cout << "  Total internal steps = " << GM.nsteps << endl;
      cout << "  Total Newton iterations = " << GM.nnewt << endl;
    }
    cout << endl;
  }



  //------ run third test: theta = 0.45 ------
  
  double theta = 0.45;
  // create time stepper object
  GeneralizedMidpoint GM(rhs, Jac, theta, y0);

  // update Newton solver parameters
  GM.newt.tol = 1.e-3;
  GM.newt.maxit = 20;
  GM.newt.show_iterates = false;
  // loop over lambda values
  for (int il=0; il<lambdas.n_elem; il++) {

    // set current lambda value into rhs and Jac objects
    rhs.lambda = lambdas(il);
    Jac.lambda = lambdas(il);

    // loop over time step sizes
    for (int ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "\nRunning generalized midpoint with stepsize h = " << h(ih)
     << ", lambda = " << lambdas(il) << ":\n";
      mat Y = GM.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      for (size_t i=0; i<Nout; i++)
        printf("    y(%.1f) = %9.6f   |error| = %.2e\n",
               tspan(i), Y(0,i), Yerr(0,i));
      Y.save("Y3.txt", raw_ascii);
      cout << "  Max error = " << Yerr.max() << endl;
      cout << "  Total internal steps = " << GM.nsteps << endl;
      cout << "  Total Newton iterations = " << GM.nnewt << endl;
    }
    cout << endl;
  }

*/


}



 