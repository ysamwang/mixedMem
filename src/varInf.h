#ifndef VARINF_H
#define VARINF_H

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "mm_model.h"
#include "settings.h"
#include "mStep.h"
#include "eStep.h"
#include "utils.h"


using namespace Rcpp;
using namespace arma;

double varInfC(mm_model model_old);
double compute_logf(mm_model model);
double compute_ELBO(mm_model model);
double alpha_Objective(mm_model model, vec alph);
double alpha_Objective(mm_model model);


#endif

