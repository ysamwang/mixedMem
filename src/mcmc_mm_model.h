#ifndef MCMC_MM_MODEL_H
#define MCMC_MM_MODEL_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "settings.h"

using namespace Rcpp ;

class mm_model_mcmc
{
public:
    //constructor
    mm_model_mcmc(List model);

    // Observed Values
    int getT();
    int getJ();
    int getK();
    int getR(int j);
    int getV(int j);
    int getObs(int i, int j, int r);
    std::string getDist(int j);

    // Get Latents
    double getLambda(int i, int k) ;
    int getZ(int i, int j, int r) ;
    //Set Latents
    void setLambda(int i, NumericVector target) ;
    void setZ(int i, int j, int r, int target) ;


    // Get Sampled Parameters
    double getTheta(int j, int k, int v);
    double getAlpha();
    double getKsi(int k);
    NumericVector getKsi();
     // Set Sampled Parameters
    void setAlpha(double target);
    void setKsi(int k, double target);
    void setKsi(NumericVector target);
    void setTheta(int j, int k, int v, double target) ;
    void setTheta(int j, int k, NumericVector target);


    // Get Hyperparameters
    double getTau(int j, int k , int v);
    double getBeta();
    double getGamma();

    // Max V for ragged arrays
    int getMaxV();
    int getMaxR();

    // Get extended model params
    int getStayerObs(int s, int j, int r);
    int getStayerMatch(int i) ;
    int getStayerStatus(int i) ;
    int getS();
    int getGomMembers();
    int getExtended();
    double getP(int s);
    // Set extended model params
    void setP(int s, double target) ;
    void setP(NumericVector target) ;
    void setStayerStatus(int i, int target) ;
    


protected:

    //Observed Quantities
    int T;  // Number of individuals
    int J;  // Number of Variables
    IntegerVector Rj ;  // Number of repeated measurements for each variable
    int K ;  // Number of sub-populations
    IntegerVector Vj ;  // Number of choices per variable
    int maxV ;  // Largest V (used for navigating ragged array)
    int maxR ;  // Largest R (used for navigating ragged array)
    CharacterVector dist ; // Distribution type of each variable
    IntegerVector obs ;  // Observed responses will be array of dimension Total x J x maxR

    //Latent Variables
    NumericVector lambda ;  // individual membership parameters (of dimension Total x K)
    NumericVector Z ; // context indicators (of dimension Total x J x maxR)

    //Parameters
    NumericVector theta ; // parameters governing distributions (of dimension J x K x maxV)
    NumericVector alpha ; // sum of membership Dirichlet parameter
    NumericVector ksi ; // relative weights of each sub-population

    //Hyper-hyper parameters
    NumericVector tau; // hyper parameters of Dirichlet/Beta prior for theta
    double beta ;  // shape parameter for alpha
    double gamma ; // scale parameter for alpha

    // Extended Objects
    int extended; // 1 for extended model; 0 for normal
    NumericVector P ; // probability of compartments. 0:(S-1) indicate stayers; S indicates GoM
    IntegerVector stayerMatch ; // Which compartment an individual could be in 0:S-1 are stayers, S is GoM
    IntegerVector stayerStatus ; // If we consider the individual a stayer or not 1 = stayer; 0 = GoM
    int S ; // Number of Stayer classes
    IntegerVector stayerObs ; // stayer signiatures of dimension (S x J x maxRj)

    // Helper functions for extended
    int checkIndStayer(int i); //
    int checkIndStayerHelp(int i, int s);

private:
};
#endif // MCMC_MM_MODEL_H

