#ifndef MM_MODELEXT_H
#define MM_MODELEXT_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "settings.h"
#include "mm_model.h"

using namespace Rcpp ;
using namespace arma;

class mm_modelExt: public mm_model
{
public:
    //constructor
    mm_modelExt(List model);
    int getFixedObs(int i, int j, int r, int n);
    double getP(int s);
    NumericVector getP();
    double getBeta(int i, int s);
    double getBeta(int s);
    NumericVector getBeta();
    int getStayersClass(int i);
    int getStayersFirstID(int s);
    NumericVector getStayers();
    double getStayersProb(int s);
    int getS();
    NumericVector getNumStayers();
    double getNumStayers(int s);

    //set individual element
    void setP(int s, double target);
    void setBeta(int s, double target);
    Rcpp::List returnModel();

protected:
    int S;
    NumericVector fixedObs;
    NumericVector stayers;
    NumericVector stayersFirstID;
    NumericVector stayersCount;
    NumericVector P;
    NumericVector beta;
    int checkIndStayer(int i);
    int checkIndStayerHelp(int i, int s);

private:
};
#endif // MM_MODEL_H
