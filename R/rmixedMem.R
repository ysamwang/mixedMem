#' Simulate Mixed Membership Data
#' 
#' Simulate data from the mixed membership generative model
#' 
#' Given a the parameters and dimensions of a mixed membership model, the function
#' returns a random sample of observed values, context indicators and group memberships.
#'    
#' @param Total the count of individuals in the sample
#' @param J the count of variables observed on each individual
#' @param Rj vector of length J specifying the number of repeated measurements
#'  on each variable
#' @param Nijr an array of dimension (Total, J, max(Rj)) indicating the number
#'  of ranking levels for each replication. For multinomial and bernoulli
#'  variables, Nijr[i,j,r] = 1. For rank variables, Nijr[i,j,r] indicates the
#'  number of candidates ranked.
#' @param K the number of sub-populations
#' @param Vj vector of length J specifying the number of possible candidates
#'  for each variable. For a bernoulli variable Vj[j] = 1.
#' @param dist vector of length J specifying variable types. Options are
#'  "bernoulli", "multinomial" or "rank" corresponing to the distributions
#'   of the observed variables
#' @param alpha vector of length K which is the hyperparameter for Dirichlet
#'  membership distribution
#' @param theta array of dimension (J,K,max(Vj)) which governs the variable
#'  distributions. theta[j,k,] is the parameter for how sub-population k responds
#'  to the variable j. If the number of candidates differs across variables, any
#'  unusued portions of theta should be 0.
#'  @param lambda a matrix containing the group membership for each individual
#' @return A list containing a random sample of lambda (group memberships),
#'  Z (context) and obs
rmixedMem <- function(Total, J, Rj, Nijr, K, Vj, dist, theta, alpha, lambda = NULL)
{
  if(is.null(lambda)){
    lambda <- gtools::rdirichlet(Total, alpha)
  }
  Z <-array(-1, dim = c(Total, J, max(Rj), max(Nijr)))
  obs <- array(-1, dim = c(Total, J, max(Rj), max(Nijr)))
  for(i in 1:Total)
  {
    for(j in 1:J)
    {
      for(r in 1:Rj[j])
      {
        if(dist[j] =="multinomial") {
          sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) # sub-population which governs the response
          obs[i,j,r,1] <- sample.int(Vj[j], size = 1, prob = theta[j,sub.pop,c(1:Vj[j])])-1 # must be in 0:(Vj[j]-1)
          Z[i,j,r,1] <- sub.pop-1
        }
        if(dist[j] == "rank") {
          select = c()
          for(n in 1:Nijr[i,j,r])
          {
            sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) # sub-population which governs the response
            prob = theta[j,sub.pop,c(1:Vj[j])]
            prob[c(select)] <- 0
            select <- c(select, sample.int(Vj[j], size = 1, replace = F, prob = prob))
            obs[i,j,r,n] <- select[n]-1  # must be in 0:(Vj[j]-1)
            Z[i,j,r,n] <- sub.pop-1
          }
        }
        if(dist[j] == "bernoulli"){
          sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) # sub-population which governs the response
          obs[i,j,r,1] <- rbinom(n = 1, size = 1, prob = theta[j,sub.pop,1])
          Z[i,j,r,1] <- sub.pop-1
        }
      }
    }
  }
  
  out <- list(lambda = lambda, Z = Z, obs = obs)
  return(out)
}