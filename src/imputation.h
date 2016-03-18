#ifndef IMPUTATION_H
#define IMPUTATION_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "mcmc_mm_model.h"
#include "utilities.h"
#include "settings.h"

void imputeLatent(mm_model_mcmc model) ;
void imputeStayers(mm_model_mcmc model) ;

double evalProb(mm_model_mcmc model, int i, int j, int k, int r) ;
NumericVector estimateGoMProb(mm_model_mcmc model, int numSamples) ;
int estimateGoMProbIndividual(mm_model_mcmc model) ;
double evalProb_condZ(mm_model_mcmc model, int i) ;

#endif // IMPUTATION_H

