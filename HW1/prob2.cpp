// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 1 -- Sep 04
// Problem 2 -- Newton's Method
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement:
// Write the Lagrange interpolating polynomial for the data
// (1,-2/3), (0,-1/6), (-1,-17/3), (2,-7/6): 
// P = x^3 - 3 * x^2 + 3/2 * x - 1/6
// Approximate the roots of this polynomial using a numerical root-finding method.
#include <iostream>
#include <cmath>

using namespace std;

// compute the value of the polynomial
double f (double x)
{   
	return pow(x,3) - 3 * pow(x,2) + (3.0/2.0) * x - 1.0/6.0;
}

// compute the value of the first derivative of the polynomial
double df (double x)
{   
    return 3 * pow(x,2) - 6 * x + 3.0/2.0;
}

int main()
{   
	double x = 0.0; // iterative x
    cout << "Enter the initial guess (recommend 0): ";
    cin >> x;
    double tmp = 0.0; // tmp = x - f(x)/f'(x)
    double *r = new double[3]; // r = roots
    double d; // difference between x(n+1) and x(n)
    for (int i = 0; i < 3; i++)
    {
    	do // Newton's Method
    	{   
    		tmp = x - (f(x)/df(x));
    		d = tmp - x;
    	    x = tmp;
    	}while(abs(f(x)) > 1e-10);
    	r[i] = x;
    	x = x + i + 1; // in order to find DIFFERENT roots
    	printf("Roots[%i] is %.10e.\n", i, r[i]);
    }
    delete [] r;
    return 0;
}
