// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 3 -- Oct 09
// Problem 3 -- Simple Implicit Methods
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
  cout << "\nRunning generalized midpoint with theta = 1.0:" << endl;
  // loop over lambda values
  for (int il=0; il<lambdas.n_elem; il++) {

    // set current lambda value into rhs and Jac objects
    rhs.lambda = lambdas(il);
    Jac.lambda = lambdas(il);

    // storage for errors at each step size
    vec e(h);
    e.fill(0.0);

    // loop over time step sizes
    for (int ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "lambda =  " << lambdas(il) << ", stepsize h = " << h(ih) << ":\n";
      mat Y = GM.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih) = Yerr.max();
      cout << "  Max error = " << Yerr.max() << endl;
    }

    mat A(h.n_elem,2);
    vec b(h.n_elem);
    for (size_t ih=0; ih<h.n_elem; ih++) {
      A(ih,0) = log(h(ih));
      A(ih,1) = 1.0;
      b(ih)   = log(e(ih));
    }
    mat AtA = A.t()*A;
    vec Atb = A.t()*b;
    vec s = solve(AtA, Atb);
    cout << "Overall convergence rate estimate = " << s(0) << endl;
    cout << endl;
  }



  //------ run second test: theta = 0.5 ------
  double theta2 = 0.5;
  // create time stepper object
  GeneralizedMidpoint GM2(rhs, Jac, theta2, y0);

  // update Newton solver parameters
  GM2.newt.tol = 1.e-3;
  GM2.newt.maxit = 20;
  GM2.newt.show_iterates = false;
  cout << "\nRunning generalized midpoint with theta = 0.5:" << endl;
  // loop over lambda values
  for (int il=0; il<lambdas.n_elem; il++) {

    // set current lambda value into rhs and Jac objects
    rhs.lambda = lambdas(il);
    Jac.lambda = lambdas(il);

    // storage for errors at each step size
    vec e(h);
    e.fill(0.0);

    // loop over time step sizes
    for (int ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "lambda =  " << lambdas(il) << ", stepsize h = " << h(ih) << ":\n";
      mat Y = GM2.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih) = Yerr.max();
      cout << "  Max error = " << Yerr.max() << endl;
    }

    mat A(h.n_elem,2);
    vec b(h.n_elem);
    for (size_t ih=0; ih<h.n_elem; ih++) {
      A(ih,0) = log(h(ih));
      A(ih,1) = 1.0;
      b(ih)   = log(e(ih));
    }
    mat AtA = A.t()*A;
    vec Atb = A.t()*b;
    vec s = solve(AtA, Atb);
    cout << "Overall convergence rate estimate = " << s(0) << endl;
    cout << endl;
  }



  //------ run third test: theta = 0.45 ------
  double theta3 = 0.45;
  // create time stepper object
  GeneralizedMidpoint GM3(rhs, Jac, theta3, y0);

  // update Newton solver parameters
  GM3.newt.tol = 1.e-3;
  GM3.newt.maxit = 20;
  GM3.newt.show_iterates = false;
  cout << "\nRunning generalized midpoint with theta = 0.45:" << endl;
  // loop over lambda values
  for (int il=0; il<lambdas.n_elem; il++) {

    // set current lambda value into rhs and Jac objects
    rhs.lambda = lambdas(il);
    Jac.lambda = lambdas(il);

    // storage for errors at each step size
    vec e(h);
    e.fill(0.0);

    // loop over time step sizes
    for (int ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "lambda =  " << lambdas(il) << ", stepsize h = " << h(ih) << ":\n";
      mat Y = GM3.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih) = Yerr.max();
      cout << "  Max error = " << Yerr.max() << endl;
    }

    mat A(h.n_elem,2);
    vec b(h.n_elem);
    for (size_t ih=0; ih<h.n_elem; ih++) {
      A(ih,0) = log(h(ih));
      A(ih,1) = 1.0;
      b(ih)   = log(e(ih));
    }
    mat AtA = A.t()*A;
    vec Atb = A.t()*b;
    vec s = solve(AtA, Atb);
    cout << "Overall convergence rate estimate = " << s(0) << endl;
    cout << endl;
  }


  printf("It is not always more accurate when the step size gets smaller. "
    "For the first method theta = 1.0, when the stiffness parameter lambda gets 10 times bigger, "
    "the max error gets 10 times smaller for each theta. It is asymptotically stable as the "
    "stability region is less than 1. \n"
    "For the second method theta = 0.5, the max errors for the same step size are almost the same "
    "with different lambda. The values of the stability region is close to 1. \n"
    "For the third method theta = 0.45, the stability region is sometimes greater than 1, "
    "which leads to unstable status. \n"
    "Again, it reminds me that the question of stability is orthogonal to the question of accuracy.\n");
  
  return 0;
}



 