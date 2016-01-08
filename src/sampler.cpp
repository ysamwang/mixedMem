#include "sampler.h"


NumericVector rDirichlet(NumericVector alpha)
{
    int i;
    double sum = 0.0;
    int v = alpha.size();
    NumericVector ret(v);
    for(i = 0 ; i < v; i ++)
    {
        ret[i] = R::rgamma(alpha[i], 1.0);
        sum += ret[i];
    }

    //calculated on log scale for numerical reasons
    for(i = 0 ; i < v; i ++)
    {
        ret[i] = exp(log(ret[i]) - log(sum));
    }
    return ret;
}



NumericVector estimateGoMProb(mm_modelExt& model, int numSamples){
    // counts of each stayer class; 0 represents stayer class
    // obsCounts(0) is the number of GoM observations
    IntegerVector obsCounts(model.getS(), 0);
    NumericVector ret(model.getS());
    int i,k;

    for(k = 0; k < numSamples; k ++)
    {
        for(i = 0; i < model.getT(); i++)
        {
            obsCounts[estimateGoMProbIndividual(model, i)]++;
        }
    }

    int s;
    int total = sum(obsCounts);

    for(s = 0; s < model.getS(); s++)
    {
        ret[s] = obsCounts[s] / total;
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


    for(j = 0; j < model.getJ(); j++)
    {


        for(r = 0; r < model.getR(j); r++)
        {

            if(model.getDist(j) == BERNOULLI)
            {
                Z = sample_cat(lambda);

                //Sample observation
                X = (int) runif(1)[0] < model.getTheta(j, Z, 0);

                //check if observation matches groups
                for(s = 1; s  < model.getS(); s++)
                {
                    if(X != model.getFixedObs(s, j, r, 0))
                    {
                        possibleClasses[s] = 0;
                    }
                }

                //Check if there are any remaining stayer groups left
                if(max(possibleClasses) == 0)
                {
                    return 0;
                }

            }
            else if(model.getDist(j) == MULTINOMIAL)
            {

                Z = sample_cat(lambda);

                //extract relevant probability vector
                NumericVector prob(model.getV(j));
                for(v = 0; v < model.getV(j); v++)
                {
                    prob[v] = model.getTheta(j, Z, v);
                }

                //Sample observation
                X = sample_cat(prob);

                //check if observation matches groups
                for(s = 1; s  < model.getS(); s++)
                {
                    if(X != model.getFixedObs(s, j, r, 0))
                    {
                        possibleClasses[s] = 0;
                    }
                }

                if(max(possibleClasses) == 0)
                {
                    return 0;
                }
            }
            else if(model.getDist(j) == RANK)
            {

            //possible items to select
             IntegerVector rankResponse(model.getV(j), 1);


                for(s = 0; s < S; s++)
                {
                    if(model.getN(i,j,r) != model.getNStayer(s,j,r))
                    {
                        possibleClasses[s] = 0;
                    }
                }

                if(max(possibleClasses) == 0)
                {
                    return 0;
                }


                for(n = 0 ; n < model.getN(i,j,r); n++)
                {
                    //extract relevant probabilities
                    Z = sample_cat(lambda);
                    int V = model.getV(j);
                    NumericVector prob(V, 0.0);
                    for(v = 0; v < V; v++)
                    {
                        if(rankResponse[v]){
                            prob[v] = model.getTheta(j,Z, v);
                        }
                    }
                    X = sample_cat(prob);

                    //X cannot be selected again later
                    rankResponse[X] = 0;

                    for(s = 0; s < S; s++)
                    {
                        if(X != model.getFixedObs(s,j,r,n))
                        {
                            possibleClasses[s] = 0;
                        }
                    }

                    if(max(possibleClasses) == 0)
                    {
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
    NumericVector p = prob / sum(prob);
    double temp = runif(1)[0];
    int ret = 0;
    while(temp > 0)
    {
        temp -= p[ret];
        ret ++;
    }
    return ret - 1;
}



