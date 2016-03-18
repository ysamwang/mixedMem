#include "mcmc_mm_model.h"

using namespace Rcpp ;

mm_model_mcmc::mm_model_mcmc(List model)
{
// Observed Values
    T = (int) as<IntegerVector>(model[0])[0];
    J = (int) as<IntegerVector>(model[1])[0];
    dist = as<CharacterVector>(model[2]);
    Rj = as<IntegerVector>(model[3]);
    Vj = as<IntegerVector>(model[4]);
    maxV = max(Vj);
    maxR = max(Rj);

    obs = as<IntegerVector>(model[5]);

// Fixed Parameter
    K = (int) as<IntegerVector>(model[6])[0];

// Sampled Parameters
    theta = clone(as<NumericVector>(model[7])) ;
    alpha = clone(as<NumericVector>(model[8])) ;
    ksi = clone(as<NumericVector>(model[9])) ;
    lambda = clone(as<NumericVector>(model[10])) ;
    Z = clone(as<NumericVector>(model[11])) ;

// Hyperparameters
    tau = as<NumericVector>(model[12]) ;
    beta = (double) as<NumericVector>(model[13])[0] ;
    gamma = (double) as<NumericVector>(model[14])[0] ;

// Extended model Parameters
    extended = (int) as<IntegerVector>(model[15])[0] ;
    P = clone(as<NumericVector>(model[16])) ;
    S = (int) as<IntegerVector>(model[17])[0] ;
    stayerObs = as<IntegerVector>(model[18]) ;
// 1 indicates in Stayer compartment, 0 indicates in GoM
    stayerStatus = IntegerVector(T, 0);

// Which stayer signiature an individual matches
// 0:(S-1) represent stayer classes
// S indicates that the individual does not match any of the stayer signiatures
    stayerMatch = IntegerVector(T, 0);
    int i;
    for(i = 0; i < T; i++) {
        stayerMatch[i] = checkIndStayer(i);
    }
}



/*
* Access model elements
*/
int mm_model_mcmc::getT()
{
    return T ;
}

int mm_model_mcmc::getJ()
{
    return J ;
}

int mm_model_mcmc::getK()
{
    return K;
}

int mm_model_mcmc::getR(int j)
{
    return Rj[j] ;
}

int mm_model_mcmc::getV(int j)
{
    return Vj[j];
}

//Objects which are arrays in R are stored a vectors in C++
// so we have get functions to help with indexing

int mm_model_mcmc::getObs(int i, int j, int r)
{
    return obs[i + T * j + (T * J)*r];
}

std::string mm_model_mcmc::getDist(int j)
{
    return as<std::string>(dist[j]);
}

double mm_model_mcmc::getLambda(int i, int k)
{
    return(lambda[i + T * k]);
}

int mm_model_mcmc::getZ(int i, int j, int r)
{
    return(Z[i + T * j + T * J * r]);
}

double mm_model_mcmc::getTheta(int j, int k, int v)
{
    return theta[j + J * k + (J * K) * v] ;
}

// Alpha is stored as an NumericVector
// so we access the first (and only) element
double mm_model_mcmc::getAlpha()
{
    return alpha[0];
}

double mm_model_mcmc::getKsi(int k)
{
    return ksi[k] ;
}

NumericVector mm_model_mcmc::getKsi()
{
    return ksi;
}

/*
* Hyper parameters
*/
double mm_model_mcmc::getTau(int j, int k , int v)
{
    return tau[j + J * k + J * K * v];
}

double mm_model_mcmc::getBeta()
{
    return beta;
}

double mm_model_mcmc::getGamma()
{
    return gamma;
}


/*
* max elements to help with ragged arrays
*/
int mm_model_mcmc::getMaxV()
{
    return maxV;
}

int mm_model_mcmc::getMaxR()
{
  return maxR;
}


/*
* Set Functions
*/
void mm_model_mcmc::setLambda(int i, NumericVector target)
{
    int k;
    for(k = 0; k < K; k++) {
        lambda[i + T * k] = target[k];
    }
}


void mm_model_mcmc::setZ(int i, int j, int r, int target)
{
    Z[i + (T * j) + (T * J * r)] = target;
}



void mm_model_mcmc::setTheta(int j, int k, int v, double target)
{
    theta[j + (J * k) + (J * K * v)] = target;
}


void mm_model_mcmc::setAlpha(double target)
{
    alpha[0] = target;
}

void mm_model_mcmc::setKsi(int k, double target)
{
    ksi[k] = target;
}

void mm_model_mcmc::setKsi(NumericVector target)
{
    int k;
    for(k = 0; k < K; k++) {
        ksi(k) = target(k);
    }
}

void mm_model_mcmc::setTheta(int j, int k, NumericVector target)
{
    int v;
    for (v = 0; v < Vj[j]; v++) {
        theta[j + J * k + J * K * v] = target[v];
    }
}

/*
* Extended GoM methods
*/
//Helper function to check whether an individual matches a stayer class
int mm_model_mcmc::checkIndStayer(int i)
{
    int s, temp;
    for(s = 0; s < S; s ++) {
        temp = checkIndStayerHelp(i, s);
        if(temp != S) {
            return temp;
        }
    }
    return S;
}

//Helper function to check whether an individual matches a stayer class
int mm_model_mcmc::checkIndStayerHelp(int i, int s)
{
    int j, r;
    for(j = 0; j < J; j++) {
        for(r = 0; r < getR(j); r++) {
            if(getObs(i, j, r) !=  getStayerObs(s, j, r)) {
                return S;
            }
        }
    }

    return s;
}

//Get Stayer Signiature
int mm_model_mcmc::getStayerObs(int s, int j, int r)
{
    return(stayerObs[s + S * j + S * J * r ]);
}

double mm_model_mcmc::getP(int s)
{
    return P[s];
}

void mm_model_mcmc::setP(int s, double target)
{
    P[s] = target;
}

void mm_model_mcmc::setP(NumericVector target)
{
    int s;
    for(s = 0; s < S + 1; s++) {
        P[s] = target[s];
    }
}

int mm_model_mcmc::getS()
{
    return S;
}

int mm_model_mcmc::getStayerStatus(int i)
{
    return stayerStatus[i];
}

void mm_model_mcmc::setStayerStatus(int i, int target)
{
    stayerStatus[i] = target;
}

int mm_model_mcmc::getStayerMatch(int i)
{
    return stayerMatch[i];
}

int mm_model_mcmc::getGomMembers()
{
    return T - sum(stayerStatus);
}

int mm_model_mcmc::getExtended()
{
        return extended;
}
