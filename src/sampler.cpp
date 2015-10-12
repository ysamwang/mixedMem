#include "sampler.h"



NumericVector rDirichlet(NumericVector alpha){
  int i;
  double sum = 0.0;
  int v = alpha.size();
  NumericVector ret(v);
  for(i = 0 ; i < v; i ++){
    ret[i] = R::rgamma(alpha[i], 1.0);
    sum += ret[i];
  }

  //calculated on log scale for numerical reasons
  for(i = 0 ; i < v; i ++){
    ret[i] = exp(log(ret[i]) - log(sum));
  }
  return ret;
}

NumericVector rMixedMem(mm_modelExt model) {
   int i,j,r;
   int Z;
   NumericMatrix lambda(model.getT(), model.getK());
   NumericVector obs;
   for(i = 0; i < model.getT(); i++){
    lambda(i,_) = rDirichlet(model.getAlpha());
     for(j = 0; j < model.getJ(); j++){
       for(r =0; r < model.getR(j); r++){
        for(n = 0; n < model.getN(i,j,r); n++){
            Z = RcppArmaillo::sample(model.getK(), 1, true, lambda(i,_));
            obs[i + model.getT() *]
            }
        }
     }
   }
 }

