#include "varInf.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List varInfInputC(Rcpp::List model_r, int print,
                    int printMod, int stepType, int maxTotalIter,
                    int maxEIter, int maxAlphaIter, int maxThetaIter, int maxLSIter,
                    double elboTol, double alphaTol, double thetaTol, double aNaught,
                    double tau, int bMax, double bNaught, double bMult, int vCutoff, SEXP holdConstSEXP) {

    mm_model model = mm_model(model_r);
    NumericVector holdConst(holdConstSEXP);
    varInfC(model, print , printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter, maxThetaIter, maxLSIter,
                              elboTol, alphaTol, thetaTol, aNaught, tau, bMax, bNaught, bMult, vCutoff, holdConst);
    return model.returnModel();
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double computeElboInputC(Rcpp::List model_r)
{
    mm_model model = mm_model(model_r);
    return compute_ELBO(model);
}


#include "mmMCMCFit.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
void mcmcInputC(Rcpp::List model_r, int burnIn, int samples, int thin, int print, Rcpp::List fileNames, int newFiles,
                  double omega, double eta, NumericVector whichWrite)
{
    mm_model_mcmc model = mm_model_mcmc(model_r);
    mixedMemMcmcC(model, burnIn, samples, thin, print, fileNames, newFiles, omega, eta, whichWrite);
}
