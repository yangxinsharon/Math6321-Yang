// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 2 -- Sep 18
// Problem 1 -- forward Euler Method
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement: Compute the true solution to the problem: 
// y'(t) = -exp(-t)*y(t), y(0) = 1.
// a) Using the forward Euler method provided in class, approximate the solution 
// to this problem numerically with stepsizes h = {0.2,0.1,0.05,0.025,0.0125} 
// for t belongs to [0,5]. Output the solution and absolute error every 1 unit of 
// time, as well as the overall maximum of these errors for the time interval.

// b) Store each of these maximum errors in a vector of the same size as your
// vector of h values.

// c) Using your stored error values, numerically estimate the order of convergence
// for Euler's method on this problem. 


#include <iostream>
#include "fwd_euler.hpp"

using namespace std;
using namespace arma;

// ODE RHS function class -- instantiates a RHSFunction
class RHS: public RHSFunction
{
	public:
		int Evaluate(double t, vec& y, vec& f)
		{    
			// evaluates the RHS function, f(t,y)
			f(0) = -exp(-t) * y(0);
    		return 0;
  		}
};

// Convenience function for analytical solution
vec ytrue(const double t)
{
	vec yt(1);
	yt(0) = exp(exp(-t)-1);
	return yt;
};


// main routine
int main()
{
	// stepsizes
	vec h("0.2 0.1 0.05 0.025 0.0125");

	// RHS initial value
	vec y0("1.0");

	// interval for time values
	double t0 = 0.0;
	double Tf = 5.0;

	// set desired output times
	int Nout = 6; // includes initial condition
	vec tspan = linspace(t0, Tf, Nout);

	// create ODE RHS function objects
	RHS f;

	// create true solution results
	mat Ytrue(1,Nout);
	for (size_t i = 0; i < Nout; i++)
	{
		Ytrue.col(i) = ytrue(tspan(i));
	}

	cout << "\nProblem 1:\n";
	ForwardEulerStepper FE(f, y0);
	vec Maxerr(size(h)); // store maximum errors in a vector of the same size of h
	for (size_t ih = 0; ih < h.n_elem; ih++)
	{
		// call stepper; output solution and error
		cout << "  h = " << h(ih) << ":\n";
		mat Y = FE.Evolve(tspan, h(ih), y0);

		// output solution, errors, and overall error
		mat Yerr = abs(Y-Ytrue);
		for (size_t i = 0; i < Nout; i++)
			printf("    y(%.1f) = %9.6f   |error| = %.2e\n",
        		tspan(i), Y(0,i), Yerr(0,i));
		cout << "  Max error = " << Yerr.max() << endl;
		Maxerr(ih) = Yerr.max();

	}

	// estimate "pointwise" order of convergence
	vec p(size(h)); // pointwise p
	double P; // overall P
	for (size_t ih = 0; ih < h.n_elem-1; ih++)
	{
		p(ih) = (log(Maxerr(ih+1))-log(Maxerr(ih))) / (log(h(ih+1))-log(h(ih)));
		cout << "  Pointwise estimate of p(" << ih << ") is " << p(ih) << endl;
		// estimate "overall" order of convergence by finding the slope of best fit line
		P = sum( (log(h(ih))-log(mean(h))) * (log(Maxerr(ih))-log(mean(Maxerr))) ) / sum( (log(h(ih))-log(mean(h))) * (log(h(ih))-log(mean(h))));
	}
	cout << "  Overall estimate of p is " << P << endl;
	return 0;
}

