#ifndef SAMPLE_ALPHA_H
#define SAMPLE_ALPHA_H

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "mcmc_mm_model.h"
#include "utilities.h"
#include "settings.h"

using namespace Rcpp;

void sampleAlpha(mm_model_mcmc model, double omega);
double propRatioAlpha(mm_model_mcmc model, double proposal, double omega);
void sampleKsi(mm_model_mcmc model, double eta);
double propRatioKsi(mm_model_mcmc model, NumericVector proposal, double eta);


#endif //SAMPLE_ALPHA_H

