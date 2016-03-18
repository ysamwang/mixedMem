#include "sampleParams.h"


void sampleTheta(mm_model_mcmc model)
{
    int i,j,r,k,v;
    int J = model.getJ();
    int K = model.getK();
    int T = model.getT();
    NumericVector target;

    for(j = 0; j < J; j ++)
    {
        if(model.getDist(j)== MULTINOMIAL)
        {

            NumericVector thetaDraw(model.getV(j));

            // Draw theta_{j,k} for each k
            for(k = 0; k < K; k++) {


                for(v = 0; v < model.getV(j); v++) {
                    thetaDraw(v) = model.getTau(j, k, v);
                }

//                Rcout << "Theta Pre: " <<thetaDraw <<std::endl;
                for(i = 0 ; i < T; i++) {
                    for(r  = 0; r < model.getR(j); r++) {
                        // Check if individual
                        if((model.getZ(i, j, r) == k) && ! model.getStayerStatus(i)) {
                            thetaDraw[model.getObs(i,j,r)] += 1.0;
                        }
                    }
                }
//            Rcout << "Theta Post: " <<thetaDraw <<std::endl;
            target = rDirichlet(thetaDraw);
            model.setTheta(j, k, target);
            }
        }
        else if (model.getDist(j)== BERNOULLI)
        {
            NumericVector gammaParam(2);
            for(k = 0; k < K; k++)
            {
                gammaParam[0] = model.getTau(j, k, 0);
                gammaParam[1] = model.getTau(j, k, 1);

//                Rcout << "Gamma Pre: " << gammaParam <<std::endl;
                for(i = 0; i < T; i++)
                {
                    for(r = 0; r < model.getR(j); r++)
                    {
                        if((model.getZ(i,j,r) == k) && !model.getStayerStatus(i))
                        {
                            //note that indexes are swapped due to parameterization of bernoulli
                                gammaParam(1 - model.getObs(i,j,r)) += 1.0;
                        }
                    }
                }
                //draw Theta
//                Rcout << "Gamma Post: " << gammaParam <<std::endl;
                double target = rbeta(1, gammaParam(0), gammaParam(1))(0);
                model.setTheta(j, k, 0, target);
            }
        }
    }
}

