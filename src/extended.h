#ifndef EXTENDED_H
#define EXTENDED_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "mm_modelExt.h"
#include "settings.h"
#include "varInfExt.h"

using namespace Rcpp;
using namespace arma;

void updateP(mm_modelExt model);
void updateBeta(mm_modelExt model);
double getStayersProb(mm_modelExt model, int s);
double getStayer_logf(mm_modelExt model, int stayerID);
void updateExt(mm_modelExt model);
//evaluates the ELBO for an individual who is a stayer
double getStayersProbI(mm_modelExt model, int i);
  

#endif

