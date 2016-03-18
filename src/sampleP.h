#ifndef SAMPLEP_H
#define SAMPLEP_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "mcmc_mm_model.h"
#include "utilities.h"
#include "settings.h"

using namespace Rcpp ;

void sampleP(mm_model_mcmc model);

#endif // SAMPLEPARAMS_H


