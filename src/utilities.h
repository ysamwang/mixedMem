#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdlib.h>
#include <math.h>
#include <Rcpp.h>


using namespace Rcpp;

NumericVector rDirichlet(NumericVector alpha);
int sampleCat_sw(NumericVector prob);

#endif //SAMPLE_ALPHA_H

