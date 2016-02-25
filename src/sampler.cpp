#include "sampler.h"


NumericVector rDirichlet(NumericVector alpha)
{
    int k;
    int K = alpha.size();
    double ret_sum = 0.0;
    NumericVector ret(K);
    for(k = 0; k < K; k++) {
        ret(k) = R::rgamma(alpha(k), 1.0);
        ret_sum += ret(k);
    }
    for(k = 0; k < K; k++) {
        ret(k) = ret(k) / ret_sum;
    }
    return ret;
}


NumericVector estimateGoMProb(mm_modelExt& model, int numSamples)
{
    // counts of each stayer class; 0 represents stayer class
    // obsCounts(0) is the number of GoM observations
    IntegerVector obsCounts(model.getS(), 0);
    NumericVector ret(model.getS());
    int i,k;

    for(k = 0; k < numSamples; k ++) {
        for(i = 0; i < model.getT(); i++) {
            int temp = estimateGoMProbIndividual(model, i);
            obsCounts[temp]++;
//            Rcout << "Class: " << temp <<std::endl;
        }
    }

    int s;
    double total = sum(obsCounts);

    for(s = 0; s < model.getS(); s++) {
        ret[s] = (double) obsCounts[s] / total;

    }
    return ret;
}

int estimateGoMProbIndividual(mm_modelExt& model, int i)
{

    //Sample membership from Alpha
    NumericVector lambda = rDirichlet(model.getAlpha());

    int j, r, n, v, s;
    int X, Z;
    int S = model.getS();

    IntegerVector frame;
    IntegerVector subPops = seq_len(model.getK()) - 1;
    IntegerVector possibleClasses = seq_len(model.getS()) - 1;


    for(j = 0; j < model.getJ(); j++) {


        for(r = 0; r < model.getR(j); r++) {

            if(model.getDist(j) == BERNOULLI) {
                Z = sample_cat(lambda);
                //Sample observation
                double temp = runif(1)[0];
                X = (temp < model.getTheta(j, Z, 0));

                //check if observation matches groups
                for(s = 1; s  < S; s++) {
                    if(X != model.getFixedObs(s, j, r, 0)) {
                        possibleClasses[s] = 0;
                    }
                }

                //Check if there are any remaining stayer groups left
                if(max(possibleClasses) == 0) {
                    return 0;
                }

            } else if(model.getDist(j) == MULTINOMIAL) {

                Z = sample_cat(lambda);

                //extract relevant probability vector
                NumericVector prob(model.getV(j));
                for(v = 0; v < model.getV(j); v++) {
                    prob[v] = model.getTheta(j, Z, v);
                }

                //Sample observation
                X = sample_cat(prob);

                //check if observation matches groups
                for(s = 1; s  < model.getS(); s++) {
                    if(X != model.getFixedObs(s, j, r, 0)) {
                        possibleClasses[s] = 0;
                    }
                }

                if(max(possibleClasses) == 0) {
                    return 0;
                }
            } else if(model.getDist(j) == RANK) {

                //possible items to select
                IntegerVector rankResponse(model.getV(j), 1);


                for(s = 0; s < S; s++) {
                    if(model.getN(i,j,r) != model.getNStayer(s,j,r)) {
                        possibleClasses[s] = 0;
                    }
                }

                if(max(possibleClasses) == 0) {
                    return 0;
                }


                for(n = 0 ; n < model.getN(i,j,r); n++) {
                    //extract relevant probabilities
                    Z = sample_cat(lambda);
                    int V = model.getV(j);
                    NumericVector prob(V, 0.0);
                    for(v = 0; v < V; v++) {
                        if(rankResponse[v]) {
                            prob[v] = model.getTheta(j,Z, v);
                        }
                    }
                    X = sample_cat(prob);

                    //X cannot be selected again later
                    rankResponse[X] = 0;

                    for(s = 0; s < S; s++) {
                        if(X != model.getFixedObs(s,j,r,n)) {
                            possibleClasses[s] = 0;
                        }
                    }

                    if(max(possibleClasses) == 0) {
                        return 0;
                    }
                }
            }
        }
    }
    return max(possibleClasses);
}


int sample_cat(NumericVector prob)
{
    double u = runif(1)[0];
    int s = 0;
    u -= prob[s];
    while(u >= 0) {
        u -= prob[++s];
    }
    return s;
}



