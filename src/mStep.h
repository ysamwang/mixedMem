#ifndef MSTEP_H
#define MSTEP_H


#define BOOST_DISABLE_ASSERTS true

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "mm_model.h"
#include "settings.h"
#include "varInf.h"



using namespace Rcpp ;
using namespace arma;
double mStep_C(mm_model model, double elbo_T, int stepType, int maxAlphaIter, int maxLSIter,
                              double alphaTol, double aNaught, double tau, NumericVector iterReached);
vec getGrad(mm_model model);
mat getHess(mm_model model);
#endif
