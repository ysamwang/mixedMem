#include "extended.h"

//Updates P based on current model estimates
void updateExt(mm_modelExt model, int method)
{
    double old_Elbo, new_Elbo;
    double tol = 10e-4;
    double conv_Crit = 1.0;
    new_Elbo = compute_ELBOExt(model);
    int count = 0;
    int s;
    // Estimates of P(X_s | alpha, theta, N)
    NumericVector sampledStayerProbs(model.getS());

    // method = 1 for ELBO
    // method = 1 for sampled
    if(method == 1){
            for(s = 0; s < model.getS(); s ++){
                sampledStayerProbs[s] = exp(getStayersProb(model, s));
            }
    } else {
        sampledStayerProbs = estimateGoMProb(model, 500);
    }

    while(conv_Crit > tol) {
        old_Elbo = new_Elbo;

        updateBeta(model, sampledStayerProbs);
        updateP(model);
        new_Elbo = compute_ELBOExt(model);
        conv_Crit = fabs((new_Elbo - old_Elbo)/old_Elbo);
        count++;
    }
}


//Updates P based on current model estimates
void updateP(mm_modelExt model)
{
    double target, total;
    total = 0.0;
    int s;
    //calculate p for all stayer classes
    for(s = 1; s < model.getS(); s++) {
        target =  model.getNumStayers(s) * model.getBeta(s) / model.getT();
        total += target;
        model.setP(s, target);
    }
    //calculate p_1 since all p's must sum to 1
    model.setP(0, 1.0 - total);
}



//Updates Beta based on current model estimates
void updateBeta(mm_modelExt model, NumericVector sampledStayerProbs)
{
    double target;
    int s;

    for(s = 1; s < model.getS(); s++) {
        target = model.getP(s) / (model.getP(0) * sampledStayerProbs[s] + model.getP(s));
        if((1.0 - target) < BUMP) {
            target = 1.0 - BUMP;
        }
        model.setBeta(s, target);
    }
}




//evaluates the ELBO for an individual who is a stayer
double getStayersProb(mm_modelExt model, int s)
{

    int stayerID = model.getStayersFirstID(s);

    double t1,t2,t3,t4;
    double phi_sum = 0.0;
    double elbo;
    int j,k,r,n;
    int K = model.getK();
    int J = model.getJ();
    double back_term;
    double dg_phi_sum;
    double phi_ik, delta_ijrnk;

    //Calculate first line and second line
    t1 = 0.0;
    t2 = 0.0;
    t3 = 0.0;
    t4 = 0.0;


    t1 = lgamma(sum(model.getAlpha())) - sum(lgamma(model.getAlpha()));

    for(k = 0; k < K; k++) {
        phi_sum += model.getPhi(stayerID,k);
    }
    dg_phi_sum = boost::math::digamma(phi_sum);

    t4 += lgamma(phi_sum);
    for(k = 0; k < K; k++) {
        phi_ik = model.getPhi(stayerID,k);
        back_term = (boost::math::digamma(phi_ik) - dg_phi_sum);
        t1 += (model.getAlpha(k) - 1.0) * back_term;

        t4 += -lgamma(phi_ik);
        t4 += (phi_ik - 1.0) * back_term;

        for(j = 0; j < J; j++) {
            for(r = 0; r < model.getR(j); r++) {
                for(n = 0; n < model.getN(stayerID, j, r); n++) {
                    delta_ijrnk = model.getDelta(stayerID, j, r, n, k);
                    t2 += delta_ijrnk*back_term;
                    t4 += delta_ijrnk*log(delta_ijrnk);
                }
            }
        }
    }

    //compute 3rd line
    t3 = getStayer_logf(model, stayerID);

    elbo = t1 + t2 + t3;// - t4;
    // Rcout <<"Elbo: " <<elbo << " t4: " <<t4 <<endl;
    return elbo;
}

//evaluates log_f of stayer (helper for getStayerProb()
double getStayer_logf(mm_modelExt model, int stayerID)
{
    double logf = 0.0;
    double back_term;
    int j, k, r, n, v, Nijr;

    for(j = 0; j < model.getJ(); j++) {
        if(model.getDist(j) == BERNOULLI) {
            n = 0;
            v = 0;
            for(r = 0; r < model.getR(j); r++) {
                for(k = 0; k < model.getK(); k++) {
                    logf += ( model.getObs(stayerID,j,r,n) ? model.getDelta(stayerID,j,r,n,k)*log(model.getTheta(j,k,v)) :
                              model.getDelta(stayerID,j,r,n,k)*log(1.0 - model.getTheta(j,k,v))) ;
                }
            }
        } //end bernoulli
        else if(model.getDist(j)== MULTINOMIAL) {
            n = 0;
            for(r = 0; r < model.getR(j); r++) {
                for(k = 0; k < model.getK(); k++) {
                    logf += model.getDelta(stayerID,j,r,n,k)*log(model.getTheta(j,k,model.getObs(stayerID,j,r,n)));
                }
            }
        } //end Multinomial
        else if(model.getDist(j) == RANK) {
            for(r = 0; r < model.getR(j); r++) {
                Nijr = model.getN(stayerID, j, r);
                for(k = 0; k < model.getK(); k++) {
                    back_term = 0.0;
                    for(n = 0; n < Nijr; n++) {
                        logf += -model.getDelta(stayerID,j,r,n,k)*log(1.0 - back_term);
                        logf += model.getDelta(stayerID,j,r,n,k)*log(model.getTheta(j,k,model.getObs(stayerID,j,r,n))) ;
                        back_term += model.getTheta(j,k,model.getObs(stayerID,j,r,n));
                    }
                }
            }
        } //end Rank
    }

    return logf;
}
