#include "mm_modelExt.h"

using namespace Rcpp ;
using namespace arma;




int mm_modelExt::getFixedObs(int i, int j, int r, int n)
{
    return(fixedObs[i + j + J*r + J*maxR*n]);
}

double mm_modelExt::getP()
{
    return P;
}

double mm_modelExt::getBeta()
{
    return beta;
}

NumericVector mm_modelExt::getStayers()
{
    return stayers;
}

int mm_modelExt::getNumStayers()
{
    return numStayers;
}


//sets list of stayers
void mm_modelExt::updateStayer() {
    int i, check;
    check = 1;
    for(i = 0; i < T; i++) {
        stayers[i] = checkIndStayer(i);
        Rcout << "Ind "<< i << " Status: " << stayers[i] <<std::endl;
//        if(check && stayers[i]){
//            stayerID = i;
//            check = 0;
//            Rcout <<"First Stayer Found!! " << i <<std::endl;
//        }
    }
//    numStayers = (int) std::accumulate(stayers.begin(), stayers.end(), 0.0);
}

// checks if an individual is a stayer
int mm_modelExt::checkIndStayer(int i)
{
    int j, r, n;
    for(j = 0; j < J; j++) {
        for(r = 0; r < getR(j); r++) {
            for(n = 0; n < getN(i, j, r); n++ ) {
                if(getObs(i, j, r, n) !=  getFixedObs(0, j, r, n)){
                    return 0;
                }
            }
        }
    }
    return 1;
}

int mm_modelExt::getStayerID()
{
    return stayerID;
}

int mm_modelExt::getStayers(int i)
{
    return stayers[i] ;
}

void mm_modelExt::setP(double target)
{
    P = target;
}

void mm_modelExt::setBeta(double target)
{
    beta = target;
}
