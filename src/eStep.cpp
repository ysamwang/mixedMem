#include "eStep.h"

using namespace boost::math;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
double eStep_C(mm_model model, double elbo_E, int maxEIter, int maxBetaIter, maxLSIter, double elboTol, double betaTol,
               double aNaught, double tau, NumericVector holdConst, NumericVector iterReached)
{
    int nE = 0;
    double old_elbo_E;
    double converged_E = 1.0; //flag for convergence


    int reset = 0;

    if(reset) {
        elbo_E = resetLocal(model);
    }


    //update local parameters and betaBar
    while ((converged_E > elboTol) && (nE < maxEIter)) {
        //UPDATE!
        //update elbo and increment E step count
        old_elbo_E = elbo_E;
        nE++;

        //update delta variational parameters
        updateDelta(model);

        //update phi variational parameters
        updatePhi(model);


        //update variational parameters for theta
        updateBetaBar(model, elbo_E, maxBetaIter, betaTol, iterReached, aNaught, tau, holdConst, iterReached);


        elbo_E = compute_ELBO(model);
        //calculate convergence criteria
        converged_E = (old_elbo_E- elbo_E)/old_elbo_E;
    }

    if (nE == maxEIter) {
//        Rcout<< "Max E Steps Reached!" <<std::endl;
        iterReached[0] = 1;
    }
    return elbo_E;
}


//referred to as  dL/ddelta in derivations
double dl_ddelta(mm_model model, int i, int j, int r, int n, int k)
{
    double dl_dd = 0.0;

    //*** BERNOULLI AND MULTINOMIAL NEED TO BE UPDATED ***

    if(model.getDist(j) == BERNOULLI) {
        dl_dd += (model.getObs(i,j,r,n) ? log(model.getTheta(j,k,0)) : log(1-model.getTheta(j,k,0)));
    } // end binomial
    else if(model.getDist(j) == MULTINOMIAL) {
        dl_dd +=log(model.getTheta(j,k,model.getObs(i,j,r,n)));
    } //end multinomial
    else if(model.getDist(j) == RANK) {
        int eta;
        double back_term = 0.0;
        dl_dd += digamma(model.getBetaBar(j,k,model.getObs(i,j,r,n)));
        //term original from the denominator that gets subtracted off
        for(eta = n; eta < model.getN(i,j,r); eta++) {
            back_term += model.getBetaBar(j,k,model.getObs(i,j,r,eta));
        }

        dl_dd -= digamma(back_term);
    }//end rank
    return dl_dd;
}



//Reset all local parameters to uniform
//option to turn on or off
double resetLocal(mm_model model)
{
    int i, j, k;
    for(i = 0; i < model.getT(); i++) {
        for(k = 0; k < model.getK(); k++ ) {
            model.setPhi(i, k, 1.0/K);
            for(j = 0; j < J; j++) {
                for(r = 0; r < model.getR(j); r++) {
                    for(n = 0; n < model.getN(i,j,r); n++) {
                        model.setDelta(i, j, r, n, k, 1.0/K);
                    }
                }
            }
        }
    }
    return compute_ELBO(model);
}

/*
* Update Delta
*/
void updateDelta(mm_model model)
{
    int i, j, k;
    int T = model.getT();
    int K = model.getK();
    double phi_sum; //helper which will be used to store sum of phi
    double phi_sum_dg;
    double delta_sum;
    double placeholder;

    for(i = 0; i < T; i++) {
        //get sum_k phi_{ik} which will be used in the delta updates
        phi_sum = 0.0;
        for(k = 0; k < K; k++) {
            phi_sum += model.getPhi(i,k);
        }

        phi_sum_dg = digamma(phi_sum); //digamma of the sum of phi's

        for(j = 0; j < J; j++) {
            for(r = 0; r < model.getR(j); r++) {
                Nijr = model.getN(i,j,r);
                for(n = 0; n < Nijr; n++) {
                    delta_sum = 0.0;
                    //update deltas
                    for(k = 0; k <  K; k++) {
                        placeholder = exp(digamma(model.getPhi(i,k)) - phi_sum_dg + dl_ddelta(model, i, j, r, n, k));
                        model.setDelta(i,j,r,n,k, placeholder);
                        delta_sum += model.getDelta(i,j,r,n,k);
                    }
                    model.normalizeDelta(i,j,r,n,delta_sum);
                }
            }
        }
    }//end update delta
}


void updatePhi(mm_model model)
{
    int i, j, n, k;
    int T = model.getT();
    int J = model.getJ();
    int K = model.getK();

    for(i = 0; i < T; i++) {
        for(k = 0; k < K; k++) {
            model.setPhi(i,k, model.getAlpha(k));
        }
    }

    for(i = 0; i < T; i++) {
        for(j = 0; j< J; j++) {
            for(r = 0; r < model.getR(j); r++) {
                Nijr = model.getN(i,j,r);
                for(n = 0; n < Nijr; n++) {
                    for(k = 0; k < K; k ++) {
                        model.incPhi(i,k, model.getDelta(i,j,r,n,k)) ;
                    }
                }
            }
        }
    }
}


double updateBetaBar(mm_model model, double elbo_E, int maxBetaIter, double betaTol, NumericVector iterReached)
{
    updateBetaBarPL(model, maxBetaIter, maxLSIter, thetaTol, aNaught, tau, holdConst, iterReached);

}

void updateBetaBarPL(mm_model model, int maxBetaIter,
                     int maxLSIter, double betaTol, double aNaught,
                     double tau, NumericVector holdConst, NumericVector iterReached)
{

    double conv_crit_m = 1.0;
    double old_obj;
    double new_obj;
    int nLS;


    /*
    * <=== Fit Beta ===>
    */

    int nB = 0;
    double a;
    for(j = 0; j < model.getJ(); j++) {
        for(k = 0; k < model.getK(); k++) {
            if(is_true( all(k != holdConst) ) ) {
                int Vj = model.getV(j);
                vec newBetaBar(Vj);
                while( (conv_crit_m > alphaTol) & (nA <= maxAlphaIter) ) {

                    old_obj = beta_Obj(model);
                    nB++;
                    //create gradient
                    grad = getGradPL(model, j,k);
                    //create hessian
                    hess = getHess(model, j, k);

                    //Note we use the negative gradient, so all updates are +!
                    update = solve(hess, -grad);

                    for(v = 0; v < Vj; v++) {
                        newBetaBar(v) = model.getBetaBar(j, k, v);
                    }
                    a = aNaught;
                    new_obj= beta_Obj(model, (newBetaBar + a * update), j, k );

                    //if the steps is already too big such that a theta is negative (which leads the objective function to be NAN)
                    //then just set ob_old to the new objective +1
                    new_obj = (any( (newBetaBar + a * update) < 0 ) ? old_obj - 1.0: new_obj);

                    nLS = 0;
                    while( (old_obj - new_obj > 0) & (nLS < maxLSIter)) {
                        nLS++;
                        a *= tau; //scale back a

                        //if in a feasible space, then check objective function
                        if(all((newBetaBar + a * update)>0)) {
                            new_obj = beta_Obj(model, (newBetaBar + a * update), j, k);
                        }
                    }

                    //update alpha in model
                    if(nLS < maxLSIter) {
                        for(v = 0; v < Vj; v++) {
                            model.setBetaBar(v, a*update(v));
                        }
                    } else {
                        new_obj = old_obj;
                    }

                    //check for convergence
                    conv_crit_m = fabs((old_obj - new_obj)/old_obj);
                }

                if( nA == maxBetaIter ) {
//            Rcout << "Max Beta Steps Reached!!" << std::endl;
                    iterReached[1] = 1;
                }
            }
        }
    } // End betaBar Update
}



vec getGradPL(mm_model model, int j, int k)
{
    int i,r,n, eta;
    int K = model.getK();
    int T = model.getT();
    int Vj = model.getV(j);
    int Nijr;
    double sum_betaBar = 0. 0;
    double tg_sum_betaBar;
    vec grad = vec(Vj);
    double t2 = 0.0;
    double diff;


    for(v = 0; v < Vj; v++) {
        diff = (model.getBeta(j, k, v) - model.getBetaBar(j, k, v));
        grad(v) = diff * polygamma(1, model.getBetaBar(j,k,v));
        t2 += diff;
    }

    double temp = t2 * polygamma(1,model.getBetaBarSum(j,k));
    for(v = 0 ; v < Vj; v++) {
        grad(v) -= temp;
    }

    for (i = 0; i < T; i ++) {
        for(r = 0; r < model.getR(j); r ++) {
            Nijr = model.getN(i, j, r);
            for(n = 0; n < Nijr; n++) {
                sum_betaBar = 0.0;
                for(eta = n; eta < Nijr; eta ++) {
                    sum_betaBar += model.getBetaBar(j, k, model.getObs(i, j, r, eta));
                }
                tg_sum_betaBar = polygamma(1,sum_betaBar) * model.getDelta(i, j, r, n, k);
                for(eta = n; eta < Nijr; eta ++) {
                    grad(model.getObs(i, j, r, eta)) -=  tg_sum_betaBar;
                }

                grad(model.getObs(i, j, r, n)) += model.getDelta(i, j, r, n, k) * polygamma(1,model.getBetaBar(j, k, model.getObs(i, j, r, n) ));

            }
        }
    }

}
return grad;
}

mat getHessPL(mm_model model, int j, int k)
{
    int i, j, r, n, k, v1, v2, v;
    int eta, eta1;
    double temp;
    double diff = 0.0;
    int Vj = model.getV(j);

    for(v = 0; v < Vj; v++) {
        diff += model.getBeta(j, k, v) - model.getBetaBar(j,k,v);
    }

    mat hess(Vj, Vj, fill::ones);
    hess = hess * (-diff * polygamma(2,model.getBetaBarSum(j, k) + polygamma(1, model.getBetaBarSum(j,k)));

    for(v = 0; v < Vj; v++) {
    hess(v,v) += (model.getBeta(j,k,v) - model.getBetaBar(j,k,v)) * polygamma(2, model.getBetaBar(j,k,v)) - polygamma(1, model.getBetaBar(j,k,v));
    }

    for(i = 0; i < model.getT(); i ++) {
    for(r = 0; r < model.getRj(j); r++) {
            for(n = 0; n < model.getN(i,j,r); n ++) {
                hess(model.getObs(i,j,r,n), model.geObs(i,j,r, n)) += model.getDelta(i, j, r, n, k) * polygamma(2, model.getBetaBar(j,k, model.getObs(i,j,r,n)));
                backEnd = 0.0;
                for(eta = n; eta < model.getN(i,j,r); eta++) {

                    backEnd+= model.getBetaBar(j,k, model.getObs(i, j, r, eta));
                }
                tetraGamma_temp = model.getDelta(i,j,r,n,k)  * polygamma(2, temp);


                for(eta = n; eta < model.getN(i, j, r); eta ++ ) {
                    hess(model.getObs(i,j,r,eta), model.geObs(i,j,r, eta)) -= tetraGamma_temp;
                    for(eta1 = n; eta1 < eta; eta1++) {
                        hess(model.getObs(i,j,r,eta), model.geObs(i,j,r, eta1)) -= tetraGamma_temp;
                        hess(model.getObs(i,j,r,eta1), model.getObs(i,j,r,eta)) -= tetraGamma_temp;
                    }

                }
            }
        }
    }

    return hess;
}



double beta_Obj(mm_model model, vec betaBar , int j, int k)
{
    double objective = 0.0;
    int i, r, Nijr, n, eta, v;
    double back_term;
    for(v = 0; v < model.getV(j); v++) {
        objective += (model.getBeta(j, k, v) - betaBar(v)) *
                     (digamma(betaBar(v))- digamma(sum(betaBar)));
        objective += log(gamma(betaBar(v)));
    }
    objective -= log(gamma(betaBar(v)));

    //Dl D Delta part
    for(i = 0; i < model.getT(); i++) {
        for(r = 0; r < model.getR(j); r++) {
            Nijr = model.getN(i,j,r);
            back_term =0.0;
            for(n = 0; n < Nijr; n++) {
                for(eta = n; eta < Nijr; eta++) {
                    back_term += betaBar(model.getObs(j,k, eta));
                }
                objective += model.getDelta(i,j,r,n,k)*(digamma(betaBar(v)) - digamma(back_term));
            }
        }
    }
    return objective;
}


