#include "extended.h"


//Updates P based on current model estimates
void updateP(mm_modelExt model){
    double target, total;
    total = 0.0;
    int s;
    //calculate p for all stayer classes
    for(s = 1; s < model.getS(); s++) {
        Rcout << "Stayers: " << model.getNumStayers(s) <<endl;
        target =  model.getNumStayers(s) * model.getBeta(s) / model.getT();
        total += target;
        model.setP(s, target); 
        Rcout <<"s: " <<s <<" assign: " << target << " beta: "<< model.getBeta(s)  <<endl;
    } 
    //calculate p_1 since all p's must sum to 1
    model.setP(0, 1.0 - total); 
}


//Updates Beta based on current model estimates
void updateBeta(mm_modelExt model) {
    double target;
    int s;
    for(s = 1; s < model.getS(); s++) {
        target = model.getP(s) / (model.getP(0) * getStayersProb(model, s) + model.getP(s));
        model.setBeta(s, target);
    }
}



//evaluates the ELBO for an individual who is a stayer
double getStayersProb(mm_modelExt model, int s){
    //ADD find first stayer
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
        t1+= (model.getAlpha(k) - 1.0) * back_term;

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

    elbo = t1 + t2 + t3 - t4;
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

