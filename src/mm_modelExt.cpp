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

  int i;
  for(i = 0; i < S; i++){
    
    if(beta[i] > (1.0 - BUMP)){
      beta[i] = (1.0 - BUMP);
    }
    stayersFirstID[i] = -1;
    stayersCount[i] = 0;
  }
  
  int temp;
  for(i = 0; i < T; i++) {
    temp = checkIndStayer(i);
    stayers[i] = temp;
    if(temp > 0 ){
      stayersCount[temp] = stayersCount[temp] + 1;
   
    }
    
    if(stayersFirstID[stayers[i]] < 0)
    {
      stayersFirstID[stayers[i]] = i;
    }
  }
}


int mm_modelExt::getFixedObs(int s, int j, int r, int n)
{
  return(fixedObs[(s-1) + (S-1)*j + (S-1)*J*r + (S-1)*J*maxR*n]);
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
  int temp = (int) stayers[i];
  if(s == 0){
    if(temp == 0) {
      return 1.0;
    } else {
      return  1.0 - beta[temp];
    }
  } else {
    if(s == temp) {
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
  
  int s, temp;
  for(s = 1; s < S; s++){
    temp = checkIndStayerHelp(i, s);
    if(temp != 0){
      return temp;
    }
  }
  return 0;
}
  
int mm_modelExt::checkIndStayerHelp(int i, int s) {
  int j, r, n;
  for(j = 0; j < J; j++) {
    for(r = 0; r < getR(j); r++) {
      for(n = 0; n < getN(i, j, r); n++ ) {
         if(getObs(i, j, r, n) !=  getFixedObs(s, j, r, n)) {
            return 0;
          }
        }
      }
    }
  return s;
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
