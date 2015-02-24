#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
const int PRINT_MOD = 1;
const int PRINT = 0;
const int SUB_PRINT = 0;


const int MAX_E_ITER = 1000; //max iterations for e step
const int MAX_A_ITER = 200; //max iterations when fitting alpha
const int MAX_TOTAL_ITER =  500; //max overall EM steps to take
const int MAX_THETA_ITER =  1000; //max steps for fitting theta
const int MAX_LS_ITER = 400; //max steps for line search

const double ELBO_TOL =  1e-6; //convergence tolerance
const double ALPHA_TOL =  1e-10; //convergence tolerance for alpha
const double THETA_TOL =  1e-10; //convergence tolerance for theta

// After normalization, if any of the values hit numerical 0, then bump by this amount
const double BUMP =  1e-16;

// Variables for Line Search
const double A_NAUGHT =  1.0; //initial step size
const double TAU =  .899; //if conditions are not met, scale step size by this factor

/*If the conditioning number for a hessian that needs to be inverted is
higher than this, then use gradient ascent instead
*/
const double COND_BOUND =  1e18;
const int B_MAX = 5; //the number of iterations in interior point method
const double B_NAUGHT = 100.0; //initial approximation to hard boundary
const double B_MULT = 100.0; //value to scale approximation to hard boundary
const int V_CUTOFF = 10; //cutoff at which to perform gradient ascent instead of inverting matrix


const std::string BERNOULLI =  "bernoulli"; //value for dist when identifying distribution type
const std::string MULTINOMIAL =  "multinomial";
const std::string RANK =  "rank";

#endif
