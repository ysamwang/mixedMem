#include "mm_modelExt.h"

using namespace Rcpp ;
using namespace arma;


mm_modelExt::mm_modelExt(List model) : mm_model::mm_model(model)
{
  fixedObs = Rcpp::clone(as<NumericVector>(model[12]));
  P = Rcpp::clone(as<NumericVector>(model[13]));
  beta =  Rcpp::clone(as<NumericVector>(model[14]));
  S =  (int) as<NumericVector>(model[15])[0];
  
  stayers = NumericVector(T);
  stayersFirstID = NumericVector(S);
  stayersCount = NumericVector(S);
  Rcout <<"Testing"<<endl;
  
  int i;
  for(i = 0; i < S; i++){
    stayersFirstID = -1;
    stayersCount = 0;
  }
  int temp;
  for(i = 0; i < T; i++) {
    temp = checkIndStayer(i);
    stayers[i] = temp;
    stayersCount[temp] = stayersCount[temp] + 1;
    
    if(stayersFirstID[stayers[i]] < 0)
    {
      stayersFirstID[stayers[i]] = i;
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
  return stayersCount;
}

double mm_modelExt::getNumStayers(int s)
{
  return stayersCount[s];
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
  IntegerVector possibleClass(S);
  for(s = 0; s < S; s++){
    possibleClass[s] = s;
  }
  
  for(j = 0; j < J; j++) {
    for(r = 0; r < getR(j); r++) {
      for(n = 0; n < getN(i, j, r); n++ ) {
        for(s = 1; s < S; s++){
          Rcout <<"fixed: " <<  getFixedObs(s-1, j, r, n) << " obs: " <<getObs(i, j, r, n);
          if(getObs(i, j, r, n) !=  getFixedObs(s-1, j, r, n)) {
            possibleClass[s] = 0;
          }
        }
      }
    }
  }
  Rcout <<Rcpp::max(possibleClass) <<endl;
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
