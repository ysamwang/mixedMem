#include "varInf.h"
#include "varInfExt.h"

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
double computeElboC(Rcpp::List model_r)
{
    mm_model model = mm_model(model_r);
    return compute_ELBO(model);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List varInfInputExtC(Rcpp::List model_r, int print,
                    int printMod, int stepType, int maxTotalIter,
                    int maxEIter, int maxAlphaIter, int maxThetaIter, int maxLSIter,
                    double elboTol, double alphaTol, double thetaTol, double aNaught,
                    double tau, int bMax, double bNaught, double bMult, int vCutoff, SEXP holdConstSEXP) {

    int s, check;
    check = 1;
    mm_modelExt model = mm_modelExt(model_r);
    for(s = 1; s < model.getS(); s++){
      if(model.getNumStayers(s) == 0)
      {
        Rcout <<"Error: No Stayers found for class " <<s <<". Terminating mmVarFit!"<<endl;
        check = 0;
      }
    }
    if(check)
    {
    NumericVector holdConst(holdConstSEXP);
    varInfExtC(model, print, printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter, maxThetaIter, maxLSIter,
                              elboTol, alphaTol, thetaTol, aNaught, tau, bMax, bNaught, bMult, vCutoff, holdConst);
    }
    return model.returnModel();
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double computeElboExtC(Rcpp::List model_r)
{
  int s, check;
  check = 1;
    mm_modelExt model = mm_modelExt(model_r);
    for(s = 1; s < model.getS(); s++){
      if(model.getNumStayers(s) == 0)
      {
        Rcout <<"Error: No Stayers found for class " <<s <<". Cannot compute ELBO!"<<endl;
        check = 0;
      }
    }
    Rcout <<"Beta: " << model.getBeta(0) <<" "<<model.getBeta(1)<<endl;
    if(check){
      return compute_ELBOExt(model);
    } else {
      return 0.0;
    }
}

//[[Rcpp::export]]
SEXP rDirichlet_sw(SEXP alpha_r){
  int i;
  NumericVector alpha = as<NumericVector>(alpha_r);
  double sum = 0.0;
  int v = alpha.size();
  NumericVector ret(v);
  for(i = 0 ; i < v; i ++){
    ret[i] = R::rgamma(alpha[i], 1.0);
    sum += ret[i];
  }

  //calculated on log scale for numerical reasons
  for(i = 0 ; i < v; i ++){
    ret[i] = exp(log(ret[i]) - log(sum));
  }
  return Rcpp::wrap(ret);
}

#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
int rsampler_sw(SEXP prob_r, SEXP sample_r){
  NumericVector prob = as<NumericVector>(prob_r);
  NumericVector s = as<NumericVector>(sample_r);
  NumericVector ret = RcppArmadillo::sample(s, 1, true, prob);
  return ret[0];
}
