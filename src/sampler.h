#ifndef SAMPLER_H
#define SAMPLER_H


#define BOOST_DISABLE_ASSERTS true

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "mm_modelExt.h"
#include "settings.h"
#include "varInfExt.h"
#include "utils.h"


using namespace Rcpp ;
using namespace arma;

NumericVector rDirichlet(NumericVector alpha);


#endif

