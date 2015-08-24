#include "mm_modelExt.h"

using namespace Rcpp ;
using namespace arma;


mm_modelExt::mm_modelExt(List model) : mm_model::mm_model(model)
{
    fixedObs = as<NumericVector>(model[12]);
    P = Rcpp::clone(as<NumericVector>(model[13]));
    beta =  Rcpp::clone(as<NumericVector>(model[14]));
    S = as<NumericVector>(model[15])[0];
    NumericVector stayers(T);
    NumericVector stayerFirstID(S);
    NumericVector stayerCount(S);
    int s;
    for(s = 0; s < S; s++)
    {
      stayerFirstID[s] = -1;
      stayerCount[s] = 0;
    }
    int i;
    for(i = 0; i < T; i++) {
        stayers[i] = checkIndStayer(i);
        stayerCount[stayers[i]] = stayerCount[stayers[i]] + 1;

        if(stayerFirstID[stayers[i]] < 0)
        {
            stayerFirstID[stayers[i]] = i;
        }
    }
}


int mm_modelExt::getFixedObs(int s, int j, int r, int n)
{
    return(fixedObs[s + j + J*r + J*maxR*n]);
}

double mm_modelExt::getP(int s)
{
    return P[s];
}

NumericVector mm_modelExt::getP()
{
    return P;
}

double mm_modelExt::getBeta(int i, int s)
{
    if(s == 0){
        if(stayers[i] == 0) {
            return 1.0;
        } else {
            return  1.0 - beta[stayers[i]];
        }
    } else {
        if(s == stayers[i] ) {
            return beta[s];
        }
    }
    return 0.0;
}

NumericVector mm_modelExt::getBeta()
{
    return beta;
}

double mm_modelExt::getBeta(int s)
{
    return beta[s];
}

NumericVector mm_modelExt::getNumStayers()
{
    return stayerCount;
}


int mm_modelExt::getS()
{
    return S;
}

NumericVector mm_modelExt::getStayers()
{
    return stayers;
}

// checks if an individual is a stayer
int mm_modelExt::checkIndStayer(int i)
{
    int j, r, n, s;
    NumericVector possibleClass(S);
    for(s = 0; s < S; s++)
    {
      possibleClass[s] = s;  
    }
    for(j = 0; j < J; j++) {
        for(r = 0; r < getR(j); r++) {
            for(n = 0; n < getN(i, j, r); n++ ) {
                for(s = 1; s < S; s++){
                    if(getObs(i, j, r, n) !=  getFixedObs(s, j, r, n)) {
                        possibleClass[s] = 0;
                    }
                }
            }
        }
    }
    return Rcpp::max(possibleClass);
}


int mm_modelExt::getStayersClass(int i)
{
    return (int) stayers[i];
}

int mm_modelExt::getStayersFirstID(int s)
{
    return (int) stayersFirstID[s];
}

void mm_modelExt::setP(int s, double target)
{
    P[s] = target;
}

void mm_modelExt::setBeta(int s, double target)
{
    beta[s] = target;
}

Rcpp::List mm_modelExt::returnModel()
{
    return Rcpp::List::create(Rcpp::Named("alpha", alpha),
                              Rcpp::Named("theta", theta),
                              Rcpp::Named("phi", phi),
                              Rcpp::Named("delta", delta),
                              Rcpp::Named("P", P),
                              Rcpp::Named("beta", beta));
}
