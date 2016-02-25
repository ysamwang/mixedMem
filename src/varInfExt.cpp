#include "varInfExt.h"

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
double varInfExtC(mm_modelExt model, int print,
               int printMod, int stepType, int maxTotalIter,
               int maxEIter, int maxAlphaIter, int maxThetaIter, int maxLSIter,
               double elboTol, double alphaTol, double thetaTol, double aNaught,
               double tau, int bMax, double bNaught, double bMult, int vCutoff, NumericVector holdConst, int method)
{



    /*
    * Initializations
    */
    NumericVector iterReached(3); //Vector to keep track of whether any of the iteration limits are hit in the E-Step or M-Steps
    double converged_T = 1.0; //convergence criteria: (old_elbo_T - elbo_T)/old_elbo_T
    double old_elbo_T = 0.0; //value of elbo from previous iteration
    double elbo_T = compute_ELBOExt(model); //updated Elbo
    int nT = 0; //count of Total EM Steps


    //only run an E-step
    if (stepType == 0) {
        NumericVector sampledStayerProbs = estimateGoMProb(model, 1000);
        updateBeta(model, sampledStayerProbs);
        elbo_T = eStepExt(model, elbo_T, maxEIter, elboTol, iterReached);
        if((nT % printMod == 0) && (print==1)) {
            Rcout<<"E-Step + Beta Update: "<<elbo_T<<endl;
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


            updateExt(model, method);
            elbo_T = compute_ELBOExt(model);
            //print if necessary
            if((nT % printMod == 0) && (print==1)) {
                Rcout<<"X-Step: "<<elbo_T<<endl;
            }



            // E Step
            elbo_T = eStepExt(model, elbo_T, maxEIter, elboTol, iterReached); //defined in eStep.cpp
            //print if necessary
            if((nT % printMod == 0) && (print==1)) {
                Rcout<<"E-Step: "<<elbo_T<<endl;
            }



            //M-step; choice of which parameters to update handled inside mStep_C function
            elbo_T = mStepExt(model, elbo_T, stepType, maxAlphaIter, maxThetaIter, maxLSIter,
                             alphaTol, thetaTol, aNaught, tau, bMax, bNaught, bMult, vCutoff, holdConst, iterReached); //defined in mStep.cpp

            //print if necessary
            if((nT % printMod == 0) && (print == 1)) {
                Rcout<<"M-Step: "<<elbo_T<<endl;
            }



            //update convergence criteria
            converged_T = fabs((old_elbo_T - elbo_T)/old_elbo_T);
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



double compute_ELBOExt(mm_modelExt model)
{
    double t1,t2,t3,t4, t5;
    double phi_sum;
    double elbo;
    int i, j, k, r, n,s;
    int T = model.getT();
    int K = model.getK();
    int J = model.getJ();
    double back_term;
    double dg_phi_sum;
    double phi_ik, delta_ijrnk;
    int Nijr;
    double fullGomMembers = T - sum(model.getBeta() * model.getNumStayers());

    //Calculate first line and second line
    t1 = 0.0;
    t2 = 0.0;
    t3 = 0.0;
    t4 = 0.0;

    t1 = fullGomMembers *(lgamma(sum(model.getAlpha())) - sum(lgamma(model.getAlpha())));

    t5 = sum( model.getNumStayers() * (model.getBeta() * log(model.getP()) + (1.0 - model.getBeta()) * log(model.getP(0))) );
    t5 += (T - sum(model.getNumStayers())) * log(model.getP(0));
    for(s = 1; s < model.getS(); s++){
      //Check at end prevents log(0) computation
      t5 += -model.getNumStayers(s) * ( model.getBeta(s) * log(model.getBeta(s)) +  (model.getBeta(s) == 1.0 ? 0.0 : ((1.0 - model.getBeta(s)) * log(1.0 - model.getBeta(s))) ) );
    }

    for(i = 0; i < T; i++) {
        phi_sum = 0.0;
        for(k = 0; k < K; k++) {
            phi_sum += model.getPhi(i, k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        t4 += lgamma(phi_sum);
        for(k = 0; k < K; k++) {
            phi_ik = model.getPhi(i, k);
            back_term = (boost::math::digamma(phi_ik) - dg_phi_sum);
            t1 += (model.getAlpha(k) - 1.0) * back_term  * model.getBeta(i, 0);

            t4 += -lgamma(phi_ik);
            t4 += (phi_ik - 1.0) * back_term;

            for(j = 0; j < J; j++) {
                for(r = 0; r < model.getR(j); r++) {
                    Nijr = model.getN(i, j, r);
                    for(n = 0; n < Nijr; n++) {
                        delta_ijrnk = model.getDelta(i, j, r, n, k);
                        t2 += delta_ijrnk * back_term * model.getBeta(i, 0);
                        t4 += delta_ijrnk * log(delta_ijrnk);
                    }
                }
            }
        }
    }

    //compute 3rd line
    t3 = compute_logfExt(model);

    elbo = t1 + t2 + t3 - t4 + t5;

    //debug!
    if(!(elbo > -INFINITY)) {
        Rcout<< t1 <<" "<<t2 <<" "<<t3 <<" "<<t4 <<" " << t5 <<endl<<"Alpha: "<<endl;
        for(k = 0; k < K; k++) {
            Rcout<<model.getAlpha(k)<<" ";
        }
        Rcout<<endl;
    }
    return(elbo);
}


double compute_logfExt(mm_modelExt model)
{
    double logf = 0.0;
    double back_term;
    int i,j,k,r,n,v, Nijr;

    for(i = 0; i < model.getT(); i++) {
        for(j = 0; j < model.getJ(); j++) {
            if(model.getDist(j) == BERNOULLI) {
                n = 0;
                v = 0;
                for(r = 0; r < model.getR(j); r++) {
                    for(k = 0; k < model.getK(); k++) {
                        logf += ( model.getObs(i,j,r,n) ? model.getDelta(i,j,r,n,k) * log(model.getTheta(j,k,v)) :
                                  model.getDelta(i,j,r,n,k)*log(1.0 - model.getTheta(j,k,v))) * model.getBeta(i, 0);
                    }
                }
            } //end bernoulli
            else if(model.getDist(j)== MULTINOMIAL) {
                n = 0;
                for(r = 0; r < model.getR(j); r++) {
                    for(k = 0; k < model.getK(); k++) {
                        logf += model.getDelta(i,j,r,n,k) * log(model.getTheta(j,k,model.getObs(i,j,r,n))) * model.getBeta(i, 0);
                    }
                }
            } //end Multinomial
            else if(model.getDist(j) == RANK) {
                for(r = 0; r < model.getR(j); r++) {
                    Nijr = model.getN(i,j,r);
                    for(k = 0; k < model.getK(); k++) {
                        back_term = 0.0;
                        for(n = 0; n < Nijr; n++) {
                            logf += -model.getDelta(i,j,r,n,k)*log(1.0 - back_term) * model.getBeta(i, 0);
                            logf += model.getDelta(i,j,r,n,k) * log(model.getTheta(j,k,model.getObs(i,j,r,n))) * model.getBeta(i, 0);
                            back_term += model.getTheta(j,k,model.getObs(i,j,r,n));
                        }
                    }
                }
            } //end Rank
        }
    }
    return(logf);
}



double alpha_ObjectiveExt(mm_modelExt model, vec alph)
{
    double objective;
    int T = model.getT();
    int K = model.getK();
    int i;
    int k;
    double back_term;
    double phi_sum;
    double dg_phi_sum;
    double sum_lgamma_alpha = 0.0;
    double fullGomMembers = T - sum(model.getBeta() * model.getNumStayers());

    for(k = 0; k < K; k++) {
        sum_lgamma_alpha += lgamma(alph(k));
    }

    objective = fullGomMembers * (lgamma(sum(alph)) - sum_lgamma_alpha);
    for(i = 0; i < T; i++) {
        phi_sum = 0.0;
        for(k = 0; k < K; k++) {
            phi_sum += model.getPhi(i,k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        for(k = 0; k < K; k++) {
            back_term = (boost::math::digamma(model.getPhi(i,k)) - dg_phi_sum);
            objective += (alph(k) - 1.0) * back_term * model.getBeta(i, 0);
        }
    }
    return objective;
}

double alpha_ObjectiveExt(mm_modelExt model)
{
    int T = model.getT();
    int K = model.getK();
    int i;
    int k;
    double back_term;
    double phi_sum;
    double dg_phi_sum;
    double fullGomMembers = T - sum(model.getBeta() * model.getNumStayers());

    double objective;
    objective = fullGomMembers * (lgamma(sum(model.getAlpha())) - sum(lgamma(model.getAlpha())));
    for(i = 0; i < T; i++) {
        phi_sum = 0.0;
        for(k = 0; k < K; k++) {
            phi_sum += model.getPhi(i,k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        for(k = 0; k < K; k++) {
            back_term = (boost::math::digamma(model.getPhi(i,k)) - dg_phi_sum);
            objective += (model.getAlpha(k) - 1.0) * back_term  * model.getBeta(i, 0);
        }
    }
    return objective;
}

