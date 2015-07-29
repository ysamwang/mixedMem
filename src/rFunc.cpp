#include "varInf.h"
#include "varInfExt.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double varInfInputC(Rcpp::List model_r, int print,
                    int printMod, int stepType, int maxTotalIter,
                    int maxEIter, int maxAlphaIter, int maxThetaIter, int maxLSIter,
                    double elboTol, double alphaTol, double thetaTol, double aNaught,
                    double tau, int bMax, double bNaught, double bMult, int vCutoff, SEXP holdConstSEXP) {

    mm_model model = mm_model(model_r);
    NumericVector holdConst(holdConstSEXP);
    double elbo = varInfC(model, print , printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter, maxThetaIter, maxLSIter,
                              elboTol, alphaTol, thetaTol, aNaught, tau, bMax, bNaught, bMult, vCutoff, holdConst);
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
double varInfInputExtC(Rcpp::List model_r, int print,
                    int printMod, int stepType, int maxTotalIter,
                    int maxEIter, int maxAlphaIter, int maxThetaIter, int maxLSIter,
                    double elboTol, double alphaTol, double thetaTol, double aNaught,
                    double tau, int bMax, double bNaught, double bMult, int vCutoff, SEXP holdConstSEXP) {

    Rcout <<" Ello Guv'nr"<<endl;
//    Rcout <<as<NumericVector>(model_r[0])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[1])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[2])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[3])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[4])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[5])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[6])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[7])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[8])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[9])[0] <<endl;
//    Rcout <<as<CharacterVector>(model_r[10])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[11])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[12])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[13])[0] <<endl;
//    Rcout <<as<NumericVector>(model_r[14])[0] <<endl;
    Rcout << "Check Done!" <<endl;

    mm_modelExt model = mm_modelExt(model_r);

    Rcout<< model.getP()<<" "<<model.getBeta()<<" "<< model.getNumStayers()<<endl;
    NumericVector holdConst(holdConstSEXP);
    Rcout << "Starting Fit"<<endl;
    double elbo = varInfExtC(model, print , printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter, maxThetaIter, maxLSIter,
                              elboTol, alphaTol, thetaTol, aNaught, tau, bMax, bNaught, bMult, vCutoff, holdConst);
    Rcout<< model.getP()<<" "<<model.getBeta()<<" "<< model.getNumStayers()<<endl;
    return(0.0);
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double computeElboExtC(Rcpp::List model_r)
{
    mm_model model = mm_modelExt(model_r);
    return compute_ELBO(model);
}
