#include "imputation.h"


using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
void imputeStayers(mm_model_mcmc model)
{
    int i, s;

//    // Calculate posterior probability of being in stayer class where P(X|GoM) is calculated by Monte Carlo
//    NumericVector sampledStayerProbs = estimateGoMProb(model, 10);
//
//    NumericVector stayerPosterior( model.getS() ) ;
//
//    for ( s = 0; s < model.getS(); s++ ) {
//        stayerPosterior[s] = model.getP(s) / (model.getP(s) + model.getP(model.getS()) * sampledStayerProbs[s] ) ;
//    }


    int status;
    // Sample stayer status from stayerPosterior
    for (i = 0; i < model.getT(); i++) {
        //If the observation matches any stayer signiature
        if ( model.getStayerMatch(i) != model.getS() ) {
            double gomProb = evalProb_condZ(model, i);
            status = runif(1)[0] < (model.getP(model.getStayerMatch(i)) /
                                    (model.getP(model.getStayerMatch(i)) + model.getP(model.getS()) * gomProb) ) ;
            model.setStayerStatus(i, status) ;
        }
    }

}


// [[Rcpp::depends(RcppArmadillo)]]
void imputeLatent(mm_model_mcmc model)
{


    // Holder for posterior parameter (multinomial) for Z
    NumericVector prob(model.getK());
    // Holder for posterior parameter (Dirichlet) of lambda
    NumericVector lambdaAlpha(model.getK());
    NumericVector newLambda(model.getK() );
    int newZ; // holder for sampled Z


    // Index and global variables
    int i, j, r, k;
    int J = model.getJ();
    int K = model.getK();
    double sum = 0.0;



    for (i = 0; i < model.getT(); i ++) {
        // Keeping track of posterior parameter for Dirichlet
        lambdaAlpha = model.getAlpha() * model.getKsi() ;

        for (j = 0; j < J; j++) {
            for(r = 0; r < model.getR(j); r ++) {

                // Compute posterior probability for each category
                sum = 0.0;
                for (k = 0; k < K; k ++) {
                    prob[k] = model.getLambda(i, k) * pow(evalProb(model, i, j, r, k), 1 - model.getStayerStatus(i));
                    sum += prob[k];
                }

                // Sample new Z based on posterior probability
                newZ = sampleCat_sw( prob / sum );

                // Update Z to new sample
                model.setZ(i, j, r, newZ) ;

                // increment lambdaAlpha the Dirichlet posterior parameter
                lambdaAlpha[ newZ ] += 1.0 ;
            }
        }

        // Sample new Lambda for individual i
        newLambda  = rDirichlet(lambdaAlpha);
        if (is_true(any(newLambda == 0 ))) {
            Rcout <<"Bad Lambda sample: " << lambdaAlpha <<std::endl;
        }
        model.setLambda(i, newLambda ) ;
    }
}

double evalProb(mm_model_mcmc model, int i, int j, int r, int k)
{
    if(model.getDist(j) == "multinomial") {
        return model.getTheta(j, k, model.getObs(i, j, r ));
    } else {

        return (model.getObs(i, j, r) ? model.getTheta(j, k, 0) : (1.0 - model.getTheta(j, k, 0)) ) ;
    }
}


double evalProb_condZ(mm_model_mcmc model, int i)
{
    int j, r;
    double ret = 1.0;
    int stayerClass = model.getStayerMatch(i);
    if(stayerClass != model.getS()) {
        for(j = 0; j < model.getJ(); j++) {
            if(model.getDist(j) == "multinomial") {
                for(r = 0; r < model.getR(j); r++) {
                    ret *= model.getTheta(j, model.getZ(i, j, r), model.getStayerObs(stayerClass, j, r ));
                }
            }  else {
                for(r = 0; r < model.getR(j); r++) {
                        if( model.getStayerObs(stayerClass, j, r )){
                          ret *= model.getTheta(j, model.getZ(i, j, r), 0);
                        } else {
                            ret *= (1.0 - model.getTheta(j, model.getZ(i, j, r), 0) );
                        }

                }
            }
        }
    }
    return ret;

}

NumericVector estimateGoMProb(mm_model_mcmc model, int numSamples)
{
    // counts of each stayer class; 0 represents stayer class
    // obsCounts(0) is the number of GoM observations
    IntegerVector obsCounts(model.getS() + 1, 0);
    NumericVector ret(model.getS() + 1);
    int i,k;
    int temp ;
    for(k = 0; k < numSamples; k++) {
        for(i = 0; i < model.getT(); i++) {
            temp = estimateGoMProbIndividual(model);
            obsCounts[temp]++;
        }
    }


    int s;
    int total = sum(obsCounts);

    for(s = 0; s < model.getS() + 1; s++) {
        ret[s] = (double) obsCounts[s] / total;
    }
    return ret;
}


int estimateGoMProbIndividual(mm_model_mcmc model)
{
    //Sample membership from Alpha


    int j, r,  v, s;
    int X, Z;

    IntegerVector frame;
    IntegerVector possibleClasses(model.getS() + 1);
    for(s = 0; s < model.getS() + 1; s ++) {
        possibleClasses[s] = s;
    }

    NumericVector lambda = rDirichlet(model.getKsi() * model.getAlpha());
    for(j = 0; j < model.getJ(); j++) {
        for(r = 0; r < model.getR(j); r++) {
            if(model.getDist(j) == BERNOULLI) {
                Z = sampleCat_sw(lambda);

                //Sample observation
                double temp = runif(1)[0];
                X = (temp < model.getTheta(j, Z, 0));
                //check if observation matches groups
                for(s = 0; s  < model.getS(); s++) {
                    if(X != model.getStayerObs(s, j, r)) {
                        possibleClasses[s] = model.getS();
                    }
                }

                //Check if there are any remaining stayer groups left
                if(min(possibleClasses) == model.getS()) {
                    return model.getS();
                }

            } else if(model.getDist(j) == MULTINOMIAL) {

                Z = sampleCat_sw(lambda);
                //extract relevant probability vector
                NumericVector prob(model.getV(j));
                for(v = 0; v < model.getV(j); v++) {
                    prob[v] = model.getTheta(j, Z, v);
                }

                //Sample observation
                X = sampleCat_sw(prob);
                //check if observation matches groups
                for(s = 0; s  < model.getS(); s++) {
                    if(X != model.getStayerObs(s, j, r)) {
                        possibleClasses[s] = model.getS();
                    }
                }

                //Check if there are any remaining stayer groups left
                if(min(possibleClasses) == model.getS()) {
                    return model.getS();
                }
            }
        }
    }

    // Gone through all observed values; check for remaining classes
    return min(possibleClasses);
}
