/* Stability-Optimized Explicit Runge-Kutta Methods:
   3 fourth-order methods: LSRK(12,4), LSRK(13,4), LSRK(14,4) with different stages
   3 test problems: 1) basic test from Niegemann
                    
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
  double alpha = 2.0;
  int Evaluate(double t, vec& y, vec& f) {    // evaluates the RHS function, f(t,y)
    f(0) = 1.0/y(0) - y(1)*exp(t*t)/(t*t) - t;
    f(1) = 1.0/y(1) - exp(t*t) - alpha*t*exp(-t*t);
    return 0;
  }
};

// Convenience function for analytical solution
vec ytrue(const double t) {
  vec yt(2);
  yt(0) = 1/t;
  yt(1) = exp(-t*t);
  return yt;
};

// main routine
int main() {

  // time steps to try
  vec h("0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 "
    "0.009 0.008 0.007 0.006 0.005 0.004 0.003 0.002 0.001 "
    "0.0009 0.0008 0.0007 0.0006 0.0005 0.0004 0.0003 0.0002 0.0001");

  // vector for errors at each h
  vec e(h.n_elem);
  e.fill(0.0);
  // vector for convergence rate
  mat conv(3,6); 
  conv.fill(0.0);

  // initial condition and time span
  vec y0("1.0 0.3678794412");
  y0(1) = exp(-1.0); // exp(-1) to be more accurate
  double t0 = 1.0;
  double Tf = 1.4;

  // set desire output times
  int Nout = 5;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);

  // create ODE RHS function object
  RHS f;
  //vec a("2.5 2.4 2.3 2.2 2.1 2.0"); //alpha

  // create true solution results
  mat Ytrue(2,Nout);
  for (size_t i=0; i<Nout; i++)
    Ytrue.col(i) = ytrue(tspan(i));
 
  // problem 1
  cout << "\nProblem 1 Basic Test from Niegemann:\n";
  //------------LSRK12---------------
  cout << "\nLSRK12:\n";
  LSRK12Stepper LSRK12(f,y0);
  // create LSRK12 solvers
  //for (size_t ia=0; ia<a.n_elem; ia++) {
    //f.alpha = a(ia);
    //cout << "  alpha = " << f.alpha << ":";

  for (size_t ih=0; ih<h.n_elem; ih++) {

    // call stepper
    cout << "  h = " << h(ih) << ":";
    mat Y = LSRK12.Evolve(tspan, h(ih), y0);

    // output steps, errors at t_end, and convergence rates
    e(ih) = abs(Y(0,4)-Ytrue(0,4)) + abs(Y(1,4)-Ytrue(1,4));
    if (ih == 9 || ih == 18 || ih == 27) {
      conv(ih/9-1,0) = h(ih);
      conv(ih/9-1,1) = log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1));
      cout << "  Total error = " << e(ih) << ",  conv rate = " 
           << log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1)) << endl;
    } else {
      cout << "  Total error = " << e(ih) << endl;
    }
  }
  //}

  e.save("basic_err_12.txt",raw_ascii);


  //------------LSRK13---------------
  e.fill(0.0);
  cout << "\nLSRK13:\n";
  LSRK13Stepper LSRK13(f,y0);
  // create LSRK13 solvers
  for (size_t ih=0; ih<h.n_elem; ih++) {

    // call stepper
    cout << "  h = " << h(ih) << ":";
    mat Y = LSRK13.Evolve(tspan, h(ih), y0);

    // output steps, errors, and convergence rates
    e(ih) = abs(Y(0,4)-Ytrue(0,4)) + abs(Y(1,4)-Ytrue(1,4));
    if (ih == 9 || ih == 18 || ih == 27) {
      conv(ih/9-1,2) = log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1));
      cout << "  Total error = " << e(ih) << ",  conv rate = " 
           << log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1)) << endl;

    } else {
      cout << "  Total error = " << e(ih) << endl;
    }
  }
  e.save("basic_err_13.txt",raw_ascii);


  //------------LSRK14---------------
  e.fill(0.0);
  cout << "\nLSRK14:\n";
  LSRK14Stepper LSRK14(f,y0);
  // create LSRK14 solvers
  for (size_t ih=0; ih<h.n_elem; ih++) {

    // call stepper
    cout << "  h = " << h(ih) << ":";
    mat Y = LSRK14.Evolve(tspan, h(ih), y0);

    // output steps, errors, and convergence rates
    e(ih) = abs(Y(0,4)-Ytrue(0,4)) + abs(Y(1,4)-Ytrue(1,4));
    if (ih == 9 || ih == 18 || ih == 27) {
      conv(ih/9-1,3) = log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1));
      cout << "  Total error = " << e(ih) << ",  conv rate = " 
           << log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1)) << endl;

    } else {
      cout << "  Total error = " << e(ih) << endl;
    }
  }
  e.save("basic_err_14.txt",raw_ascii);


  //------------ERK4---------------
  e.fill(0.0);
  cout << "\nERK4:\n";
  ERK4Stepper ERK4(f,y0);
  // create LSRK14 solvers
  for (size_t ih=0; ih<h.n_elem; ih++) {

    // call stepper
    cout << "  h = " << h(ih) << ":";
    mat Y = ERK4.Evolve(tspan, h(ih), y0);

    // output steps, errors, and convergence rates
    e(ih) = abs(Y(0,4)-Ytrue(0,4)) + abs(Y(1,4)-Ytrue(1,4));
    if (ih == 9 || ih == 18 || ih == 27) {
      conv(ih/9-1,4) = log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1));
      cout << "  Total error = " << e(ih) << ",  conv rate = " 
           << log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1)) << endl;

    } else {
      cout << "  Total error = " << e(ih) << endl;
    }
  }
  e.save("basic_err_erk4.txt",raw_ascii);


//------------FE---------------
  e.fill(0.0);
  cout << "\nFE:\n";
  ForwardEulerStepper FE(f,y0);
  // create forward euler solvers
  for (size_t ih=0; ih<h.n_elem; ih++) {

    // call stepper
    cout << "  h = " << h(ih) << ":";
    mat Y = FE.Evolve(tspan, h(ih), y0);

    // output steps, errors, and convergence rates
    e(ih) = abs(Y(0,4)-Ytrue(0,4)) + abs(Y(1,4)-Ytrue(1,4));
    if (ih == 9 || ih == 18 || ih == 27) {
      conv(ih/9-1,5) = log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1));
      cout << "  Total error = " << e(ih) << ",  conv rate = " 
           << log(e(ih)/e(ih-1))/log(h(ih)/h(ih-1)) << endl;

    } else {
      cout << "  Total error = " << e(ih) << endl;
    }
  }
  e.save("basic_err_fe.txt",raw_ascii);



  h.save("basic_h.txt",raw_ascii);
  conv.save("basic_conv.txt",raw_ascii);

  return 0;
}
