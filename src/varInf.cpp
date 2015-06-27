
#include "varInf.h"

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
double varInfC(mm_model model, int print,
               int printMod, int stepType, int maxTotalIter,
               int maxEIter, int maxAlphaIter, int maxBetaIter, int maxLSIter,
               double elboTol, double alphaTol, double betaTol, double aNaught,
               double tau, NumericVector holdConst)
{

    /*
    * Initializations
    */
    NumericVector iterReached(3); //Vector to keep track of whether any of the iteration limits are hit in the E-Step or M-Steps
    double converged_T = 1.0; //convergence criteria: (old_elbo_T - elbo_T)/old_elbo_T
    double old_elbo_T = 0.0; //value of elbo from previous iteration
    double elbo_T = compute_ELBO(model); //updated Elbo
    int k; //indexing variable
    int nT = 0; //count of Total EM Steps

    //only run an E-step
    if (stepType == 0) {
        elbo_T = eStep_C(model, elbo_T, maxEIter, maxBetaIter, maxLSIter, elboTol, betaTol,
                         aNaught, tau, holdConst, iterReached, stepType);
        if((nT % printMod == 0) && (print==1)) {
            Rcout<<"E-Step: "<<elbo_T<<endl;
        }
    } else {     //go through E and M steps
        while((converged_T > elboTol) && (nT < maxTotalIter)) {

            nT++; //increment count of iterations

            //print if necessary
            if((nT % printMod == 0) && (print == 1)) {
                Rcout<<"Iter: "<<nT<<" Elbo: "<< elbo_T<<endl;
            }

            //store old Elbo
            old_elbo_T = elbo_T;

            // E Step
            elbo_T = eStep_C(model, elbo_T, maxEIter, maxBetaIter, maxLSIter, elboTol, betaTol,
                             aNaught, tau, holdConst, iterReached, stepType); //defined in eStep.cpp

            //print if necessary
            if((nT % printMod == 0) && (print==1)) {
                Rcout<<"E-Step: "<<elbo_T<<endl;
            }

            if(stepType == 3) {
                //M-step; choice of which parameters to update handled inside mStep_C function
                elbo_T = mStep_C(model, elbo_T, stepType, maxAlphaIter, maxLSIter,
                                 alphaTol, aNaught, tau, iterReached); //defined in mStep.cpp
                if((nT % printMod == 0) && (print == 1)) {
                    Rcout<<"M-Step: "<<elbo_T<<endl;
                    for(k =0; k < model.getK(); k++) {
                        Rcout<<model.getAlpha(k)<<" ";
                    }
                    Rcout<<endl;
                }

                //print if necessary

            }

            //update convergence criteria
            converged_T = (old_elbo_T - elbo_T)/old_elbo_T;
        }
    }

    //print results
    Rcout <<"Fit Complete! Elbo: " <<elbo_T<< " Iter: " << nT<<endl;

    //check if any of the iteration limits were hit
    if (nT == maxTotalIter) {
        Rcout<< "Warning: Max Total Iterations Reached!" <<endl;
    }

    if( iterReached[0] == 1) {
        Rcout<< "Warning: Max E-Step Iterations Reached!" <<endl;
    }

    if( iterReached[1] == 1) {
        Rcout << "Warning: Max Alpha Iterations Reached!" <<endl;
    }

    if ( iterReached[2] == 1 ) {
        Rcout<< "Warning: Max Theta Iterations Reached!" <<endl;
    }
    return(elbo_T);
}



double compute_ELBO(mm_model model)
{
    double t1, t2, t3, t4, t5;
    int i, j, k, v, n, r;
    int T = model.getT();
    int J = model.getJ();
    int K = model.getK();
    double sumAlpha = sum(model.getAlpha());
    t1 = 0.0;
    for(j = 0; j < J; j++) {
        for(k = 0; k < K; k++) {
            t1 += lgamma(model.getBetaSum(j, k)) - lgamma(model.getBetaBarSum(j, k));
            for(v = 0; v < model.getV(j); v++) {
                t1 += -lgamma(model.getBeta(j,k,v)) + lgamma(model.getBetaBar(j,k,v));
                t1 += (model.getBeta(j,k,v) - model.getBetaBar(j,k,v)) * (digamma(model.getBetaBar(j,k,v)) - digamma(model.getBetaBarSum(j,k)));
            }
        }
    }

    //Note uses overloaded lgamma; first one takes doubles, second one is vectorized sugar
    t2 = T * lgamma(sumAlpha) - T*sum(Rcpp::lgamma(model.getAlpha()));
    t3 = 0.0;
    t4 = 0.0;
    t5 = 0.0;

    double phiSum;
    double backTerm;
    double delta;
    for(i = 0; i < T; i++) {

        phiSum = 0.0;
        for(k = 0; k < K; k++) {
            phiSum += model.getPhi(i,k);
        }
        t5 += -lgamma(phiSum);
        phiSum = digamma(phiSum);

        for(k = 0; k < K; k++) {
            t5 += lgamma(model.getPhi(i,k));
            backTerm = digamma(model.getPhi(i,k)) - phiSum;
            t2 += (model.getAlpha(k) - 1.0) * backTerm;
            t5 += -(model.getPhi(i,k) - 1.0) * backTerm;
            for(j = 0; j < J; j++) {
                for(r = 0; r < model.getR(j); r++) {
                    for(n = 0; n < model.getN(i,j,r); n++) {
                        delta = model.getDelta(i,j,r,n,k);
                        t3 += delta * backTerm;
                        t5 += -delta * log(delta);
                    }
                }
            }

        }
    }

    t4 = compute_logf(model);

    double elbo = t1+ t2 + t3 + t4 + t5;
    return elbo;
}



double compute_logf(mm_model model)
{
    double logf = 0.0;
    double back_term;
    int i,j,k,r,n, Nijr;

    for(i = 0; i < model.getT(); i++) {
        for(j = 0; j < model.getJ(); j++) {
            if(model.getDist(j) == BERNOULLI) {
//                n = 0;
//                v = 0;
//                for(r = 0; r < model.getR(j); r++) {
//                    for(k = 0; k < model.getK(); k++) {
//                        logf += ( model.getObs(i,j,r,n) ? model.getDelta(i,j,r,n,k)*log(model.getTheta(j,k,v)) :
//                                  model.getDelta(i,j,r,n,k)*log(1.0 - model.getTheta(j,k,v))) ;
//                    }
//                }
            } //end bernoulli
            else if(model.getDist(j)== MULTINOMIAL) {
//                n = 0;
//                for(r = 0; r < model.getR(j); r++) {
//                    for(k = 0; k < model.getK(); k++) {
//                        logf += model.getDelta(i,j,r,n,k)*log(model.getTheta(j,k,model.getObs(i,j,r,n)));
//                    }
//                }
            } //end Multinomial
            else if(model.getDist(j) == RANK) {
                for(r = 0; r < model.getR(j); r++) {
                    Nijr = model.getN(i,j,r);
                    for(k = 0; k < model.getK(); k++) {
                        back_term = model.getBetaBarSum(j,k);
                        for(n = 0; n < Nijr; n++) {
                            logf += model.getDelta(i,j,r,n,k) *
                                    (digamma(model.getBetaBar(j, k, model.getObs(i,j,r,n))) - digamma(back_term));
                            back_term += -model.getBetaBar(j, k, model.getObs(i,j,r,n) );
                        }
                    }
                }
            } //end Rank
        }
    }
    return(logf);
}


double alpha_Objective(mm_model model, vec alph)
{
    double objective;
    int T = model.getT();
    int K = model.getK();
    int i;
    int k;
    double phi_sum;
    double dg_phi_sum;
    double sum_lgamma_alpha = 0.0;
    double alphSum = sum(alph);

    for(k = 0; k < K; k++) {
        sum_lgamma_alpha += lgamma(alph(k));
    }

    objective = T*lgamma(alphSum) - T*sum_lgamma_alpha;
    for(i = 0; i < T; i++) {
        phi_sum = 0.0;
        for(k = 0; k < K; k++) {
            phi_sum += model.getPhi(i,k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        for(k = 0; k < K; k++) {
            objective += (alph(k)-1)*(boost::math::digamma(model.getPhi(i,k)) - dg_phi_sum);
        }
    }
    return objective;
}

double alpha_Objective(mm_model model)
{
    int T = model.getT();
    int K = model.getK();
    int i;
    int k;
    double phi_sum;
    double dg_phi_sum;
    double alphaSum = sum(model.getAlpha());

    double objective;
    objective = T*lgamma(alphaSum) - T*sum(Rcpp::lgamma(model.getAlpha()));
    for(i = 0; i < T; i++) {
        phi_sum = 0.0;
        for(k = 0; k < K; k++) {
            phi_sum += model.getPhi(i,k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        for(k = 0; k < K; k++) {
            objective += (model.getAlpha(k)-1)*(boost::math::digamma(model.getPhi(i,k)) - dg_phi_sum);
        }
    }
    return objective;
}
