// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 1 -- Sep 04
// Problem 3 -- Newton's Method
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement:
// Write a program that uses Newton's method to find a root to the
// nonlinear system of equations x^2 + y^2 = 4; 1 = xy;
// using a relative tolerance of 1e-8, an absolute tolerance of 1e-12, 
// and with an initial guess of (x, y) = (1, 2).
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

//*****************************************************************************//
int main ()
{
	vec tmp(2); // tmp = x(n+1) - x(n)
	mat Jac(2,2,fill::zeros); // Jacobian matrix
	vec F(2); // F  = Jac * dx
	vec X(2); // initial guess
	X(0) = 1; // x
	X(1) = 2; // y
	//double *R[4]; // R = roots 
	//R[0] = new double[2];
	//R[1] = new double[2];
	//R[2] = new double[2];
	//R[3] = new double[2];
	double *R = new double[2]; // to find a root
	double relative_tol;
	double absolute_tol;
	//for (int i = 0; i < 4; i++)
	//{
		do // Newton's Method
		{
			Jac(0,0) = 2 * X(0);
			Jac(0,1) = 2 * X(1);
			Jac(1,0) = X(1);
			Jac(1,1) = X(0);
			F(0) = (pow(X(0),2) + pow(X(1),2) - 4) * (-1);
			F(1) = (X(0) * X(1) - 1) * (-1);
			tmp = solve(Jac,F);
			X = X + tmp;
			relative_tol = sqrt(tmp(0)*tmp(0) + tmp(1)*tmp(1)) / sqrt(X(0)*X(0) + X(1)*X(1));
			absolute_tol = sqrt(tmp(0)*tmp(0) + tmp(1)*tmp(1));
		}while(relative_tol >= 1e-8 || absolute_tol >= 1e-12);
		//R[i][0] = X(0);
		//R[i][1] = X(1);
		R[0] = X(0);
		R[1] = X(1);
		//printf("Roots[%i] is (%.10e,%.10e).\n", i, R[i][0],R[i][1]);
		//printf("Roots[%i] is (%.10e,%.10e).\n", i, R[0],R[1]);
		printf("Root is (%.10e,%.10e).\n", R[0],R[1]);
	//}

	//delete [] R[0];
	//delete [] R[1];
	//delete [] R[2];
	//delete [] R[3];
	delete [] R;
	return 0;
}