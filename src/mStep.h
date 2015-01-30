#ifndef MSTEP_H
#define MSTEP_H

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "mm_model.h"
#include "settings.h"
#include "varInf.h"
#include "utils.h"



using namespace Rcpp ;
using namespace arma;

double mStep_C(mm_model model, int print, double elbo);
void updateTheta(mm_model model);
vec getGrad(mm_model model);
mat getHess(mm_model model);
vec getGradPL(mm_model model, int j, int k, double b);
mat getHessPL(mm_model model, int j, int k, double b);
double rank_Objective(mm_model model, vec theta, int j, int k, double b);
void update_PL_Theta(mm_model model, int j);

#endif
