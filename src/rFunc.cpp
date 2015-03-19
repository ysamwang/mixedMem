#include "varInf.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double varInfInputC(Rcpp::List model_r)
{
    mm_model model = mm_model(model_r);
    double elbo = varInfC(model);
    return(elbo);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double computeElboC(Rcpp::List model_r)
{
    mm_model model = mm_model(model_r);
    return compute_ELBO(model);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double heldOutInputC(Rcpp::List model_r)
{
    mm_model model = mm_model(model_r);
    double elbo = heldOutC(model);
    return(elbo);
}
