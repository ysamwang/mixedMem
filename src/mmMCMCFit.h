#ifndef MCMC_MISSING_H
#define MCMC_MISSING_H


#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <RcppArmadillo.h>
#include "mcmc_mm_model.h"
#include "imputation.h"
#include "sampleParams.h"
#include "sampleAlpha.h"
#include "sampleP.h"
#include "settings.h"

void mixedMemMcmcC(mm_model_mcmc model, int burnIn, int sample, int thin, int print,
                   Rcpp::List fileNames, int newFiles, double omega,
                    double eta, NumericVector which_write);
void writeState(mm_model_mcmc model, Rcpp::List fileNames, NumericVector which_write);
void clearFiles(Rcpp::List fileNames, NumericVector which_write);

#endif // MCMC_MISSING_H

