#ifndef MM_MODEL_H
#define MM_MODEL_H

#define BOOST_DISABLE_ASSERTS true

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>

using namespace Rcpp ;
using namespace arma;

class mm_model
{
public:
    //constructor
    mm_model(List model);
    int indN(int i, int j, int r);
    int indPhi(int i, int k);
    int indDelta(int i, int j, int r, int n, int k);
    int indObs(int i, int j, int r, int n);
    int indBeta(int j, int k, int v);

    //get individual element
    NumericVector getAlpha();
    double getAlpha(int k);
    int getT();
    int getJ();
    int getK();
    int getR(int j);
    int getV(int j);
    int getN(int i, int j, int r);
    double getPhi(int i, int k);
    double getDelta(int i, int j, int r, int n, int k);
    double getBeta(int j, int k, int v);
    double getBetaBar(int j, int k, int v);
    double getBetaBarSum(int j, int k);
    double getBetaSum(int j, int k);
    int getObs(int i, int j, int r, int n);
    CharacterVector getDist();
    std::string getDist(int j);

    //set individual element
    void setAlpha(int k, double target);
    void setPhi(int i, int k, double target);
    void setDelta(int i, int j, int r, int n, int k, double target);
    void setBetaBar(int j, int k, int v, double target);

    //update helpers
    void normalizeDelta(int i, int j, int r, int n, double delta_sum);
    void incPhi(int i, int k, double inc);
    void incAlpha(int k, double inc);


protected:
    int T;
    int J;
    IntegerVector Rj ;
    int maxR ;
    IntegerVector Nijr ;
    int maxN ;
    int K ;
    IntegerVector Vj ;
    int maxV ;
    NumericVector alpha;
    NumericVector phi ;
    NumericVector delta ;
    NumericVector obs ;
    CharacterVector dist;
    NumericVector betaBar;
    NumericVector beta;
private:
};
#endif // MM_MODEL_H
