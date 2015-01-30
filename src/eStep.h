#ifndef ESTEP_H
#define ESTEP_H

#include<stdlib.h>
#include<math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "mm_model.h"
#include "settings.h"
#include "varInf.h"
#include "utils.h"


using namespace Rcpp ;
using namespace arma;

double eStep_C(mm_model model, int print, double elbo_E);
double dl_ddelta(mm_model model, int i, int j, int r, int n, int k);

#endif

