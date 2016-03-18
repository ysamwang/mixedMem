#ifndef SAMPLEPARAMS_H
#define SAMPLEPARAMS_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "mcmc_mm_model.h"
#include "utilities.h"
#include "settings.h"

using namespace Rcpp ;

void sampleTheta(mm_model_mcmc model);

#endif // SAMPLEPARAMS_H


