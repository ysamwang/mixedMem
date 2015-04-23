#' Compute the approximate BIC
#' 
#' Computes the approximate BIC of a given \code{mixedMemModel}, where the lower bound on the log-likelihood
#' (also called ELBO) is used instead of the intractable true log-likelihood. 
#' 
#' \eqn{BIC = -2 ELBO + p \log(Total)}
#' 
#' where p is the number of estimated parameters and Total is the number of individuals in
#' the sample. 
#'
#' This BIC model selection criteria is used
#' in Erosheva (2007). When counting the number of fitted parameters, we only include the
#' hyperparameters \eqn{\theta} and \eqn{\alpha} and omit the variational parameters \eqn{\phi} and \eqn{\delta}.
#'  
#' @param model the \code{mixedMemModel} object for which the BIC will be calculated
#' @return the approximate BIC value
#' @references 
#' Erosheva, E. A., Fienberg, S. E., & Joutard, C. (2007). Describing disability through individual-level mixture models for multivariate binary data. The annals of applied statistics, 1(2), 346.
#' @export
computeBIC= function(model)
{
  elbo = computeELBO(model)
  num_param = model$K
  for(j in 1:model$J)
  {
    if(model$dist[j] =="binomial")
    {
      num_param = num_param + model$K
    } else if(model$dist[j] == "multinomial")
    {
      num_param = num_param + model$K*(model$Vj[j]-1)
    } else if(model$dist[j] == "rank")
    {
      num_param = num_param + model$K*(model$Vj[j]-1)
    }
  }
  return(-2*elbo+num_param*log(model$Total))
}