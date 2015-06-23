#include "mStep.h"


using namespace boost::math;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

vec getGrad(mm_model model)
{
    int i,c,k;
    int K = model.getK();
    int T = model.getT();
    double sum_phi = 0.0;
    vec grad = vec(K);

    double sum_alpha = sum(model.getAlpha());
    double dg_sum_alpha = digamma(sum_alpha);
    for(k = 0; k < K; k++) {
        grad(k) = T*(dg_sum_alpha - digamma(model.getAlpha(k)));
        for(i = 0; i < T; i++) {
            sum_phi = 0.0;
            for(c = 0; c < K; c++) {
                sum_phi += model.getPhi(i,c);
            }
            grad(k) += digamma(model.getPhi(i,k)) - digamma(sum_phi);
        }
    }
    return grad;
}

mat getHess(mm_model model)
{
    int k;
    int K = model.getK();
    double tri_gam_eval;
    mat hess = mat(K,K);
    double sum_alpha = sum(model.getAlpha());
    tri_gam_eval = trigamma(sum_alpha);
    hess.ones();
    hess = hess*tri_gam_eval*model.getT();
    for(k = 0; k < K; k++) {
        hess(k,k) -= trigamma(model.getAlpha(k))*model.getT();
    }
    return hess;
}

double mStep_C(mm_model model, double elbo_T, int stepType, int maxAlphaIter, int maxLSIter,
               double alphaTol, double aNaught, double tau, NumericVector iterReached)
{
    int K = model.getK();
    vec grad = vec(K);
    mat hess = mat(K, K);
    vec update = vec(K);
    vec new_alpha = vec(K);
    int k;


    double conv_crit_m = 1.0;
    double old_obj;
    double new_obj;
    int nLS;

    /*
    * <=== Fit Alpha ===>
    */

    int nA = 0;
    double a;

    //Update Alpha
    if (stepType == 2 || stepType == 3) {
        while( (conv_crit_m > alphaTol) & (nA <= maxAlphaIter) ) {
            old_obj = alpha_Objective(model);
            nA++;

            //create gradient
            grad = getGrad(model);
            //create hessian
            hess = getHess(model);

            update = solve(hess,-grad);

            for(k = 0; k < K; k++) {
                new_alpha(k) = model.getAlpha(k);
            }
            a = aNaught;
            new_obj= alpha_Objective(model, (new_alpha + a * update) );

            //if the steps is already too big such that a theta is negative (which leads the objective function to be NAN)
            //then just set ob_old to the new objective +1
            new_obj = (any( (new_alpha + a * update) < 0 ) ? old_obj - 1.0: new_obj);

            nLS = 0;
            while( (old_obj - new_obj > 0) & (nLS < maxLSIter)) {
                nLS++;
                a *= tau; //scale back a

                //if in a feasible space, then check objective function
                if(all((new_alpha + a * update)>0)) {
                    new_obj = alpha_Objective(model, (new_alpha + a * update));
                }
            }

            //update alpha in model
            if(nLS < maxLSIter) {
                for(k = 0; k < K; k++) {
                    model.incAlpha(k, a*update(k));
                }
            } else {
                new_obj = alpha_Objective(model);
            }

            //check for convergence
            conv_crit_m = fabs((old_obj - new_obj)/old_obj);
        }

        if( nA == maxAlphaIter ) {
            Rcout << "Max Alpha Steps Reached!!" << std::endl;
            iterReached[1] = 1;
        }
    } // End Alpha Update
    double elbo = compute_ELBO(model);
    return elbo;
} //end m-step
