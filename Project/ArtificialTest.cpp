/* Stability-Optimized Explicit Runge-Kutta Methods:
   3 fourth-order methods: LSRK(12,4), LSRK(13,4), LSRK(14,4) with different stages
   3 test problems: 2) artificial test problem from Cash
                       increasing stiff as lambda increases
                    
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
  double lambda;
  int Evaluate(double t, vec& y, vec& f) {    // evaluates the RHS function, f(t,y)
    f(0) = -lambda*y(0) + (lambda-1.0)*exp(-t);
    return 0;
  }
};

// Convenience function for analytical solution
vec ytrue(const double t) {
  vec yt(1);
  yt(0) = exp(-t);
  return yt;
};

// main routine
int main() {

  // time steps to try
  vec h("0.1 0.05 0.04 0.03 0.02 0.01 0.005 0.001");

  // initial condition and time span
  vec y0("1.0");
  double t0 = 0.0;
  double Tf = 20.0;
  vec lambdas("100.0 200.0 400.0");

  // matrix for errors at each h and lambdas
  mat e(h.n_elem,lambdas.n_elem);
  e.fill(0.0);
  mat conv(h.n_elem-1,lambdas.n_elem);
  conv.fill(0.0);

  // set desired output times
  int Nout = 21;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);

  // create ODE RHS function object
  RHS f;

  // create true solution results
  mat Ytrue(1,Nout);
  for (size_t i=0; i<Nout; i++)
    Ytrue.col(i) = ytrue(tspan(i));

  // problem 2
  cout << "\nProblem 2 Artificial Problem from Cash:\n";
  //------------LSRK12---------------
  cout << "\nLSRK12:\n";
  for (size_t il=0; il<lambdas.n_elem; il++) {
    f.lambda = lambdas(il);
    cout << "  lambda = " << f.lambda << ":\n";
    LSRK12Stepper LSRK12(f,y0);
    // create LSRK12 solvers
    for (size_t ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "  h = " << h(ih) << ":";
      mat Y = LSRK12.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih,il) = Yerr.max();
      if (ih > 0) {
        conv(ih-1,il) = log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1));
        cout << "  Max error = " << e(ih,il) << ",  conv rate = " 
             << log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1)) << endl;
      } 
      else {
        cout << "  Max error = " << e(ih,il) << endl;
      } 

    }
    e.save("artificial_err_12.txt",raw_ascii);
    conv.save("artificial_conv_12.txt",raw_ascii);
  }


  
  //------------LSRK13---------------
  e.fill(0.0);
  conv.fill(0.0);
  cout << "\nLSRK13:\n";
  for (size_t il=0; il<lambdas.n_elem; il++) {
    f.lambda = lambdas(il);
    cout << "  lambda = " << f.lambda << ":\n";
    LSRK13Stepper LSRK13(f,y0);
    // create LSRK13 solvers
    for (size_t ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "  h = " << h(ih) << ":";
      mat Y = LSRK13.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih,il) = Yerr.max();
      if (ih > 0) {
        conv(ih-1,il) = log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1));
        cout << "  Max error = " << e(ih,il) << ",  conv rate = " 
             << log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1)) << endl;
      } 
      else {
        cout << "  Max error = " << e(ih,il) << endl;
      } 

    }
    e.save("artificial_err_13.txt",raw_ascii);
    conv.save("artificial_conv_13.txt",raw_ascii);
  }

 
  //------------LSRK14---------------
  e.fill(0.0);
  conv.fill(0.0);
  cout << "\nLSRK14:\n";
  for (size_t il=0; il<lambdas.n_elem; il++) {
    f.lambda = lambdas(il);
    cout << "  lambda = " << f.lambda << ":\n";
    LSRK14Stepper LSRK14(f,y0);
    // create LSRK14 solvers
    for (size_t ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "  h = " << h(ih) << ":";
      mat Y = LSRK14.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih,il) = Yerr.max();
      if (ih > 0) {
        conv(ih-1,il) = log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1));
        cout << "  Max error = " << e(ih,il) << ",  conv rate = " 
             << log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1)) << endl;
      } 
      else {
        cout << "  Max error = " << e(ih,il) << endl;
      } 

    }
    e.save("artificial_err_14.txt",raw_ascii);
    conv.save("artificial_conv_14.txt",raw_ascii);
  }

  //------------ERK4---------------
  e.fill(0.0);
  conv.fill(0.0);
  cout << "\nERK4:\n";
  for (size_t il=0; il<lambdas.n_elem; il++) {
    f.lambda = lambdas(il);
    cout << "  lambda = " << f.lambda << ":\n";
    ERK4Stepper ERK4(f,y0);
    // create ERK4 solvers
    for (size_t ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "  h = " << h(ih) << ":";
      mat Y = ERK4.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih,il) = Yerr.max();
      if (ih > 0) {
        conv(ih-1,il) = log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1));
        cout << "  Max error = " << e(ih,il) << ",  conv rate = " 
             << log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1)) << endl;
      } 
      else {
        cout << "  Max error = " << e(ih,il) << endl;
      } 

    }
    e.save("artificial_err_erk4.txt",raw_ascii);
    conv.save("artificial_conv_erk4.txt",raw_ascii);
  }
 


//------------FE---------------
  e.fill(0.0);
  conv.fill(0.0);
  cout << "\nFE:\n";
  for (size_t il=0; il<lambdas.n_elem; il++) {
    f.lambda = lambdas(il);
    cout << "  lambda = " << f.lambda << ":\n";
    ForwardEulerStepper FE(f,y0);
    // create ERK4 solvers
    for (size_t ih=0; ih<h.n_elem; ih++) {

      // call stepper
      cout << "  h = " << h(ih) << ":";
      mat Y = FE.Evolve(tspan, h(ih), y0);

      // output solution, errors, and overall error
      mat Yerr = abs(Y-Ytrue);
      e(ih,il) = Yerr.max();
      if (ih > 0) {
        conv(ih-1,il) = log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1));
        cout << "  Max error = " << e(ih,il) << ",  conv rate = " 
             << log(e(ih,il)/e(ih-1,il))/log(h(ih)/h(ih-1)) << endl;
      } 
      else {
        cout << "  Max error = " << e(ih,il) << endl;
      } 

    }
    e.save("artificial_err_fe.txt",raw_ascii);
    conv.save("artificial_conv_fe.txt",raw_ascii);
  } 

  h.save("artificial_h.txt",raw_ascii);
  
  return 0;
}
