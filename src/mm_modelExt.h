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
    mm_modelExt(List model) : mm_model(model)
    {
        fixedObs = as<NumericVector>(model[12]);
        P = (double) as<NumericVector>(model[13])[0];
        beta = (double) as<NumericVector>(model[14])[0];
        NumericVector stayers(T);
        stayerID = 0;
        updateStayer();
        Rcout << "New Constructor Complete!"<<std::endl;
    }

    int getFixedObs(int i, int j, int r, int n);
    double getP();
    double getBeta();
    NumericVector getStayers();
    int getStayers(int i);
    int getNumStayers();
    double getStayerProb();
    int getStayerID();
    void updateStayer();

    //set individual element
    void setP(double target);
    void setBeta(double target);


protected:
    NumericVector fixedObs;
    NumericVector stayers;
    double P;
    double beta;
    int numStayers;
    int stayerID;


    int checkIndStayer(int i);

private:
};
#endif // MM_MODEL_H
