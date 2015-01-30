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
const int MAX_THETA_ITER =  1000;
const int MAX_LS_ITER = 400;

const double ELBO_TOL =  1e-6; //convergence tolerance
const double ALPHA_TOL =  1e-10;
const double THETA_TOL =  1e-10;

const double BUMP =  1e-100;

const double A_NAUGHT =  1.0;
const double TAU =  .899;

const double COND_BOUND =  1e18;
const int B_MAX = 5;
const double B_NAUGHT = 100.0;
const double B_MULT = 100.0;
const int V_CUTOFF = 10;


const std::string BERNOULLI =  "bernoulli"; //value for dist when identifying distribution type
const std::string MULTINOMIAL =  "multinomial";
const std::string RANK =  "rank";

#endif
