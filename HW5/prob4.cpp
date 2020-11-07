// Southern Methodist University -- Math Department
// Math 6321 -- Fall 2020
// Homework 5 -- Nov 06
// Problem 4 -- Hamiltonian systems of ODEs
// Sharon Yang -- xiny@smu.edu


//*****************************************************************************//
// Problem Statement: Compute the deviation in the Hamiltonian using
// (a) the 'classical' RK4 method
// (b) the fifth-order Radau IA method
// (c) the 2-stage Gauss IRK method
// (d) the AdaptRKF solver from HW 4

#include <iostream>
#include "adapt_rkf.hpp"
#include "erk4.hpp"
#include "irk.hpp"

using namespace std;
using namespace arma;


//    ODE RHS function class -- instantiates a RHSFunction
class MyRHS: public RHSFunction {
public:
  int Evaluate(double t, vec& y, vec& f) {    // evaluates the RHS function, f(t,y)
    double D = 0.0378652;
    double S = 1.814;
    double q0 = 1.41;
    double M = 0.9953;
    f(0) = y(1)/M;					// Mq' = p 
    f(1) = -2.0*D * (1.0-exp(-S*(y(0)-q0))) * S * exp(-S*(y(0)-q0));		// p' = f(q) = -grad U(q)
    return 0;
  }
};


//    ODE RHS Jacobian function class -- instantiates a RHSJacobian
class MyJac: public RHSJacobian {
public:
  double D = 0.0378652;
  double S = 1.814;
  double q0 = 1.41;
  double M = 0.9953;
  int Evaluate(double t, vec& y, mat& J) {    // evaluates the RHS Jacobian, J(t,y)
    J(0,0) = 0.0;
    J(0,1) = 1/M;
    J(1,0) = -2.0*D*S * (2.0*S* pow(exp(-S*(y(0)-q0)),2) -  S* exp(-S*(y(0)-q0)));
    J(1,1) = 0.0;
    return 0;
  }
};


//    Hamiltonian function
double Ham(double q, double p, vec& info) {
	double H;
	double t1;
	double t2;
	double t3;
	double t4;
	t1 = p*p/info(3)/2.0;
	t2 = exp(-info(1)*(q-info(2)));
	t3 = 1.0-t2;
	t4 = info(0)*pow(t3,2);

	H = t1 + t4;
	return H;
}



// main routine
int main() {

  // time steps to try
  vec h("2, 0.2");

  // set problem information
  vec y0 = {1.4155, 1.545*0.9953/48.888};
  double t0 = 0.0;
  double Tf = 2000.0;
  vec info = ("0.0378652, 1.814, 1.41, 0.9953"); // store D,S,q0,M

  // set desired output times
  int Nout = 1001;  // includes initial condition
  vec tspan = linspace(t0, Tf, Nout);

  tspan.save("tspan.txt", raw_ascii);
  // create ODE RHS and Jacobian objects
  MyRHS rhs;
  MyJac Jac;

  double H0 = Ham(y0(0),y0(1),info);


  //---- ERK4 ---- 
  ERK4Stepper ERK4(rhs, y0);
  mat H_devi(Nout,2);
  for (size_t ih=0; ih<h.n_elem; ih++) {
    mat Y = ERK4.Evolve(tspan, h(ih), y0);

    for (size_t i=0; i<Nout; i++) {
		H_devi(i,ih) = Ham(Y(0,i),Y(1,i),info) - H0; 
    }
    Y.save("Y_ERK4.txt", raw_ascii);
  }

  H_devi.save("H_devi_ERK4.txt", raw_ascii);



  //////// RadauIA 3 stage method -- O(h^5) accurate ////////

  // create IRK stepper object
  mat RIA3_A(3,3);
  vec RIA3_b(3), RIA3_c(3);
  RIA3_A(0,0) = 1.0/9.0;
  RIA3_A(0,1) = -(1.0 + sqrt(6.0))/18.0;
  RIA3_A(0,2) = -(1.0 - sqrt(6.0))/18.0;
  RIA3_A(1,0) = 1.0/9.0;
  RIA3_A(1,1) = (88.0 + 7.0*sqrt(6.0))/360.0;
  RIA3_A(1,2) = (88.0 - 43.0*sqrt(6.0))/360.0;
  RIA3_A(2,0) = 1.0/9.0;
  RIA3_A(2,1) = (88.0 + 43.0*sqrt(6.0))/360.0;
  RIA3_A(2,2) = (88.0 - 7.0*sqrt(6.0))/360.0;
  RIA3_b(0) = 1.0/9.0;
  RIA3_b(1) = (16.0 + sqrt(6.0))/36.0;
  RIA3_b(2) = (16.0 - sqrt(6.0))/36.0;
  RIA3_c(0) = 0.0;
  RIA3_c(1) = (6.0 - sqrt(6.0))/10.0;
  RIA3_c(2) = (6.0 + sqrt(6.0))/10.0;
  IRKStepper RIA3(rhs, Jac, y0, RIA3_A, RIA3_b, RIA3_c);


  RIA3.newt.tol = 1e-3;
  RIA3.newt.maxit = 20;
  RIA3.newt.show_iterates = false;
  H_devi.fill(0.0);
  for (int ih=0; ih<h.n_elem; ih++) {
    // call stepper
    mat Y = RIA3.Evolve(tspan, h(ih), y0);
    for (size_t i=0; i<Nout; i++) {
    	H_devi(i,ih) = Ham(Y(0,i),Y(1,i),info) - H0; 
    }
    Y.save("Y_RIA3.txt", raw_ascii);
  }
  H_devi.save("H_devi_RIA3.txt", raw_ascii);




  //////// Gauss-Legendre 2 stage method -- O(h^4) accurate ////////

  // create IRK stepper object
  mat GL2_A(2,2);
  vec GL2_b(2), GL2_c(2);
  GL2_A(0,0) = 0.25;
  GL2_A(0,1) = (3.0-2.0*sqrt(3.0))/12.0;
  GL2_A(1,0) = (3.0+2.0*sqrt(3.0))/12.0;
  GL2_A(1,1) = 0.25;
  GL2_b.fill(0.5);
  GL2_c(0) = (3.0 - sqrt(3.0))/6.0;
  GL2_c(1) = (3.0 + sqrt(3.0))/6.0;
  IRKStepper GL2(rhs, Jac, y0, GL2_A, GL2_b, GL2_c);

  GL2.newt.tol = 1e-3;
  GL2.newt.maxit = 20;
  GL2.newt.show_iterates = false;
  H_devi.fill(0.0);
  for (int ih=0; ih<h.n_elem; ih++) {
    // call stepper
    mat Y = GL2.Evolve(tspan, h(ih), y0);
    for (size_t i=0; i<Nout; i++) {
    	H_devi(i,ih) = Ham(Y(0,i),Y(1,i),info) - H0; 
    }
    Y.save("Y_GL2.txt", raw_ascii);
  }
  H_devi.save("H_devi_GL2.txt", raw_ascii);




  //---- AdaptRKF ----
  // tolerances
  vec rtols("1.e-4, 1.e-6");
  vec atols(2);  atols.fill(1.e-12); 
  H_devi.zeros();
  AdaptRKF RKF(rhs, 0.0, atols, y0);
  for (size_t ir=0; ir<rtols.n_elem; ir++) {
  	RKF.rtol = rtols(ir);
  	mat Y = RKF.Evolve(tspan, y0);
    for (size_t i=0; i<Nout; i++) {
    	H_devi(i,ir) = Ham(Y(0,i),Y(1,i),info) - H0; 
    }
    Y.save("Y_RKF.txt", raw_ascii);
  }
  H_devi.save("H_devi_RKF.txt", raw_ascii);

  return 0;
}
