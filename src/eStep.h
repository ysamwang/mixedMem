#ifndef ESTEP_H
#define ESTEP_H


#define BOOST_DISABLE_ASSERTS true

#include<stdlib.h>
#include<math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "mm_model.h"
#include "settings.h"
#include "varInf.h"


using namespace Rcpp ;
using namespace arma;

double eStep_C(mm_model model, double elbo_E, int maxEIter, int maxBetaIter, int maxLSIter, double elboTol, double betaTol,
               double aNaught, double tau, NumericVector holdConst, NumericVector iterReached, int stepType);
void updatePhi(mm_model model);
void updateDelta(mm_model model);
void updateBetaBar(mm_model model, double elbo_E, int maxBetaIter, double betaTol,
                     int maxLSIter, double aNaught, double tau,
                     NumericVector holdConst, NumericVector iterReached);
double resetLocal(mm_model model);
double dl_ddelta(mm_model model, int i, int j, int r, int n, int k);
void updateBetaBarPL(mm_model model, int maxBetaIter,
                     int maxLSIter, double betaTol, double aNaught,
                     double tau, NumericVector holdConst, NumericVector iterReached);
vec getGradPL(mm_model model, int j, int k);
mat getHessPL(mm_model model, int j, int k);
double beta_Obj(mm_model model, vec betaBar , int j, int k);
double beta_Obj(mm_model model, int j, int k);
#endif

