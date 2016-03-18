#include "utilities.h"

NumericVector rDirichlet(NumericVector alpha)
{
    int k;
    int K = alpha.size();
    double ret_sum = 0.0;
    NumericVector ret(K);
    for(k = 0; k < K; k++) {
        // for numerical stability, we don't want values too close to 0
        ret(k) = std::max(R::rgamma(alpha(k), 1.0), 1e-10);
        ret_sum += ret(k);
    }

    for(k = 0; k < K; k++) {
        ret(k) = ret(k) / ret_sum;
    }

    return ret;
}

int sampleCat_sw(NumericVector prob)
{
    double u = runif(1)[0];
    int s = 0;
    u -= prob[s];
    while(u >= 0) {
        u -= prob[++s];
    }
    return s;
}

