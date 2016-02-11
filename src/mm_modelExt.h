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
    int getNumStayers(int s);

    //set individual element
    void setP(int s, double target);
    void setBeta(int s, double target);
    Rcpp::List returnModel();
    int getNStayer(int s, int j, int r);

protected:
    // Number of s-classes (including GoM class)
    int S;

    // Signiature of stayer classes
    NumericVector fixedObs;

    // vector that contains the s-class of each individual; 0 indicates GoM
    NumericVector stayers;

    // vector of length S - 1; contains the first occuring possible stayer for each s-class
    NumericVector stayersFirstID;

    //vector of length S; contains the number of individuals in each s-class; does not count GoM classes
    NumericVector stayersCount;

    //vector of length S; estimates for P
    NumericVector P;

    //vector of length S; posterior estimates for each stayer class
    NumericVector beta;

    // N(i,j,r) for each stayer signiature
    NumericVector NStayers;

    int checkIndStayer(int i);
    int checkIndStayerHelp(int i, int s);

private:
};
#endif // MM_MODEL_H
