#include "varInf.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double varInfInputC(Rcpp::List model_r, int print,
                    int printMod, int stepType, int maxTotalIter,
                    int maxEIter, int maxAlphaIter, int maxBetaIter, int maxLSIter,
                    double elboTol, double alphaTol, double betaTol, double aNaught,
                    double tau, SEXP holdConstSEXP) {

    mm_model model = mm_model(model_r);
    NumericVector holdConst(holdConstSEXP);
    double elbo = varInfC(model, print, printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter, maxBetaIter, maxLSIter,
                              elboTol, alphaTol, betaTol, aNaught, tau, holdConst);
    return(elbo);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double computeElboC(Rcpp::List model_r)
{
    mm_model model = mm_model(model_r);
    return compute_ELBO(model);
}

