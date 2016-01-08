#ifndef SAMPLER_H
#define SAMPLER_H


#define BOOST_DISABLE_ASSERTS true

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "mm_modelExt.h"
#include "settings.h"
//
//
using namespace Rcpp ;
using namespace arma;

NumericVector rDirichlet(NumericVector alpha);
NumericVector estimateGoMProb(mm_modelExt& model, int numSamples);
int estimateGoMProbIndividual(mm_modelExt& model, int i);
int sample_cat(NumericVector prob);

#endif

