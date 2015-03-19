#include "varInf.h"

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
double varInfC(mm_model model)
{
    /*
    * Initializations
    */
    double converged_T = 1.0;
    double old_elbo_T = 0.0;
    double elbo_T = compute_ELBO(model);
    int k,nT = 0; //count of Total EM Steps

    while((converged_T > ELBO_TOL) && (nT < MAX_TOTAL_ITER))
    {
        nT++;
        if((nT % PRINT_MOD == 0) && (PRINT ==1))
        {
            Rcout<<"Iter: "<<nT<<" Elbo: "<< elbo_T<<endl;
        }
        old_elbo_T = elbo_T;
        elbo_T = eStep_C(model, PRINT, elbo_T);
        if((nT % PRINT_MOD == 0) && (PRINT==1))
        {
            Rcout<<"E-Step: "<<elbo_T<<endl;
        }
        elbo_T = mStep_C(model, PRINT, elbo_T);
        if((nT % PRINT_MOD == 0) && (PRINT ==1))
        {
            Rcout<<"M-Step: "<<elbo_T<<endl;

            //UPDATE
            for(k =0; k < model.getK(); k++)
            {
                Rcout<<model.getAlpha(k)<<" ";
            }
            Rcout<<endl;
        }


        converged_T = (old_elbo_T- elbo_T)/old_elbo_T;
    }
    Rcout <<"Fit Complete! Elbo: " <<elbo_T<< " Iter: " << nT<<endl;
    return(elbo_T);
}

double heldOutC(mm_model model)
{
    double elbo_T;
    elbo_T = eStep_C(model, PRINT, 0.0);
    return(elbo_T);
}



double compute_ELBO(mm_model model)
{
    double t1,t2,t3,t4;
    double phi_sum;
    double elbo;
    int i,j,k,r,n;
    int T = model.getT();
    int K = model.getK();
    int J = model.getJ();
    double back_term;
    double dg_phi_sum;
    double phi_ik, delta_ijrnk;
    int Nijr;


    //Calculate first line and second line
    t1 = 0.0;
    t2 = 0.0;
    t3 = 0.0;
    t4 = 0.0;


    t1 = T*lgamma(sum(model.getAlpha())) - T*sum(lgamma(model.getAlpha()));
    for(i = 0; i < T; i++)
    {
        phi_sum = 0.0;
        for(k = 0; k < K; k++)
        {
            phi_sum += model.getPhi(i,k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        t4 += lgamma(phi_sum);
        for(k = 0; k < K; k++)
        {
            phi_ik = model.getPhi(i,k);
            back_term = (boost::math::digamma(phi_ik) - dg_phi_sum);
            t1+= (model.getAlpha(k)-1)*back_term;

            t4 += -lgamma(phi_ik);
            t4 += (phi_ik-1)*back_term;

            for(j = 0; j < J; j++)
            {
                for(r = 0; r < model.getR(j); r++)
                {
                    Nijr = model.getN(i,j,r);
                    for(n = 0; n < Nijr; n++)
                    {
                        delta_ijrnk = model.getDelta(i,j,r,n,k);
                        t2 += delta_ijrnk*back_term;
                        t4 += delta_ijrnk*log(delta_ijrnk);
                    }
                }
            }
        }
    }

    //compute 3rd line
    t3 = compute_logf(model);

    elbo = t1+t2+t3-t4;

    //debug!
    if(!(elbo > -INFINITY))
    {
        Rcout<< t1 <<" "<<t2 <<" "<<t3 <<" "<<t4 <<endl<<"Alpha: "<<endl;
        for(k = 0; k < K; k++)
        {
            Rcout<<model.getAlpha(k)<<" ";
        }
        Rcout<<endl;
    }
    return(elbo);
}


double compute_logf(mm_model model)
{
    double logf = 0.0;
    double back_term;
    int i,j,k,r,n,v, Nijr;

    for(i = 0; i < model.getT(); i++)
    {
        for(j = 0; j < model.getJ(); j++)
        {
            if(model.getDist(j) == BERNOULLI)
            {
                n = 0;
                v = 0;
                for(r = 0; r < model.getR(j); r++)
                {
                    for(k = 0; k < model.getK(); k++)
                    {
                        logf += ( model.getObs(i,j,r,n) ? model.getDelta(i,j,r,n,k)*log(model.getTheta(j,k,v)) :
                                  model.getDelta(i,j,r,n,k)*log(1.0 - model.getTheta(j,k,v))) ;
                    }
                }
            } //end bernoulli
            else if(model.getDist(j)== MULTINOMIAL)
            {
                n = 0;
                for(r = 0; r < model.getR(j); r++)
                {
                    for(k = 0; k < model.getK(); k++)
                    {
                        logf += model.getDelta(i,j,r,n,k)*log(model.getTheta(j,k,model.getObs(i,j,r,n)));
                    }
                }
            } //end Multinomial
            else if(model.getDist(j) == RANK)
            {
                for(r = 0; r < model.getR(j); r++)
                {
                    Nijr = model.getN(i,j,r);
                    for(k = 0; k < model.getK(); k++)
                    {
                        back_term = 0.0;
                        for(n = 0; n < Nijr; n++)
                        {
                            logf += -model.getDelta(i,j,r,n,k)*log(1.0 - back_term);
                            logf += model.getDelta(i,j,r,n,k)*log(model.getTheta(j,k,model.getObs(i,j,r,n))) ;
                            back_term += model.getTheta(j,k,model.getObs(i,j,r,n));
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
    double back_term;
    double phi_sum;
    double dg_phi_sum;
    double sum_lgamma_alpha = 0.0;

    for(k = 0; k < K; k++)
    {
        sum_lgamma_alpha += lgamma(alph(k));
    }

    objective = T*lgamma(sum(alph)) - T*sum_lgamma_alpha;
    for(i = 0; i < T; i++)
    {
        phi_sum = 0.0;
        for(k = 0; k < K; k++)
        {
            phi_sum += model.getPhi(i,k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        for(k = 0; k < K; k++)
        {
            back_term = (boost::math::digamma(model.getPhi(i,k)) - dg_phi_sum);
            objective += (alph(k)-1)*back_term;
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
    double back_term;
    double phi_sum;
    double dg_phi_sum;

    double objective;
    objective = T*lgamma(sum(model.getAlpha())) - T*sum(lgamma(model.getAlpha()));
    for(i = 0; i < T; i++)
    {
        phi_sum = 0.0;
        for(k = 0; k < K; k++)
        {
            phi_sum += model.getPhi(i,k);
        }
        dg_phi_sum = boost::math::digamma(phi_sum);

        for(k = 0; k < K; k++)
        {
            back_term = (boost::math::digamma(model.getPhi(i,k)) - dg_phi_sum);
            objective += (model.getAlpha(k)-1)*back_term;
        }
    }
    return objective;
}


