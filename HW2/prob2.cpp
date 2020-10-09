// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 2 -- Sep 18
// Problem 2 -- forward Euler Method
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement: y1' = alpha - y1 - (4*y1*y2)/(1+pow(y1,2));
// y2' = beta * y1 * (1 - y2/(1+pow(y1,2)));
// a) Set alpha = 10 and beta = 2. Run the forward Euler method using a step size
// of h = e-3, with initial condition y = [0;2], over the time interval [0,40].
// Store the solution every 0.1 units of time, and output the corresponding t,
// y1(t) and y2(t) values to disk.
// time, as well as the overall maximum of these errors for the time interval.

// b) Repeat (a) for beta = 4.
// c) Repeat (a) for beta = 3.55, and time interval [0,100].


#include <iostream>
#include "fwd_euler.hpp"

using namespace std;
using namespace arma;

static const double alpha = 10.0;
static const double beta = 3.55; // 2

// Problem 2
// ODE RHS function class -- instantiates a RHSFunction
class RHS2: public RHSFunction
{
public:
  int Evaluate(double t, vec& y, vec& f)
  { 
    // evaluates the RHS function, f(t,y)
    f(0) = alpha - y(0) - (4.0*y(0)*y(1) / (1.0+pow(y(0),2)));
    f(1) = beta*y(0) * (1.0 - y(1)/(1.0+pow(y(0),2)));
    return 0;
  }
};


// main routine
int main()
{
  // File for store and output results
  FILE *outFile = fopen("2c_results.txt","w");


  // time steps to try
  vec h("0.001");

  // set problem information
  vec y0_2("0.0 2.0");
  double t0 = 0.0;
  double Tf = 100.0; //40.0

  // set desired output times
  int Nout = 1001;  // includes initial condition // 401
  vec tspan = linspace(t0, Tf, Nout);

  // create ODE RHS function objects
  RHS2 f2;

  // problem 2: loop over time step sizes; call stepper and compute errors
  cout << "\nProblem 2:\n";
  ForwardEulerStepper FE2(f2, y0_2);
  for (size_t ih = 0; ih < h.n_elem; ih++)
  {
    // call stepper; output solution and error
    cout << "  h = " << h(ih) << ":\n";
    mat Y = FE2.Evolve(tspan, h(ih), y0_2); // store the solution
    fprintf(outFile, " t        y1(t)        y2(t)\n");
    for (size_t i = 0; i < Nout; i++)
    {
      printf("    y1(%.1f) = %9.6f   y2(%.1f) = %9.6f\n",
             tspan(i), Y(0,i), tspan(i), Y(1,i));
      fprintf(outFile, "%.1f    %9.6f    %9.6f\n", 
      tspan(i), Y(0,i), Y(1,i));
      //fprintf(outFile, "t = %.1f    y1(%.1f) = %9.6f    y2(%.1f) = %9.6f\n", 
      //tspan(i), tspan(i), Y(0,i), tspan(i), Y(1,i));
    }
    fclose(outFile);

  }
  return 0;
}
