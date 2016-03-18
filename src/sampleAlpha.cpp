#include "sampleAlpha.h"

void sampleAlpha(mm_model_mcmc model, double omega)
{
    double proposal;
    double alpha = model.getAlpha();
    proposal = R::rgamma(omega, alpha / omega);
    double propRatio = propRatioAlpha(model, proposal, omega);
    if (log(runif(1)[0]) < propRatio) {
        model.setAlpha(proposal);
    }
}

// Returns the log prop ratio
double propRatioAlpha(mm_model_mcmc model, double proposal, double omega)
{
    double ret = 0.0;
    int i, k;
    double log_ksi;
    double inner_sum = 0.0;
    double alpha = model.getAlpha();

    for(k = 0; k < model.getK(); k++) {
        log_ksi = 0.0;
        for(i = 0; i < model.getT(); i++) {
            log_ksi += (1.0 - model.getStayerStatus(i)) * log(model.getLambda(i,k));
            if (!std::isfinite(log(model.getLambda(i,k)))) {
//                Rcout <<"Lambda: " << model.getLambda(i,k)<< std::endl;
            }
        }
        inner_sum += model.getKsi(k) * log_ksi;
    }



    ret = (log(proposal) - log(alpha)) * (model.getBeta() - 1.0);
    if (!std::isfinite(ret)) {
        Rcout <<"Ret 1"<<std::endl;
    }

    ret += -(model.getGamma() - inner_sum) * ( proposal - alpha );
    if (!std::isfinite(ret)) {
        Rcout <<"Ret 2"<<std::endl;
    }

    ret += model.getGomMembers() * ((lgamma(proposal) - sum(lgamma( (model.getKsi() * proposal) ))) - (lgamma(alpha) - sum(lgamma((model.getKsi() * alpha)))));
    if (!std::isfinite(ret)) {
        Rcout <<"Ret 3"<<std::endl;
    }
    ret += (log(alpha) - log(proposal)) * (2.0 * omega - 1.0) - omega * (alpha / proposal - proposal / alpha);
    if (!std::isfinite(ret)) {
        Rcout <<"Ret 4"<<std::endl;
    }
    return ret;
}



void sampleKsi(mm_model_mcmc model, double eta)
{
    NumericVector propAlpha = eta * model.getK() * model.getKsi();
    NumericVector proposal = rDirichlet(propAlpha);

    double propRatio = propRatioKsi(model, proposal, eta);

    if(log(runif(1)(0)) < propRatio) {
        model.setKsi(proposal);
    }
}

double propRatioKsi(mm_model_mcmc model, NumericVector proposal, double eta)
{
    int i, k;
    double ret = 0.0;
    double alpha = model.getAlpha();
    double log_lambda = 0.0;
    int K = model.getK();

    for(k = 0; k < K; k++) {
        log_lambda = 0.0;
        for(i = 0; i < model.getT(); i++) {
            log_lambda += (1.0 - model.getStayerStatus(i)) * log(model.getLambda(i,k));
        }
        ret += log_lambda *(proposal(k) - model.getKsi(k)) * alpha;
    }

    ret += model.getGomMembers() * (sum(lgamma( model.getKsi() * alpha)) - sum((lgamma(proposal * alpha))));
    ret += sum(lgamma(eta * K * model.getKsi())) - sum(lgamma(eta * K * proposal));
    ret += sum(log(model.getKsi()) * (proposal - 1.0));
    ret += -sum(log(proposal) * (model.getKsi() - 1.0));
    return ret;
}

