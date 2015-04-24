#' Constructor for a Mixed Membership Model
#' 
#' Constructor for a \code{mixedMemModel} object which can be used for analysis 
#' in the mixedMem package.
#' 
#' The function returns an object of \code{mixedMemModel} class. This object holds dimensions of the model,
#' the observed data, and the estimates of the model parameters. Once a \code{mixedMemModel} object is created,
#' it can be fit using the \code{mmVarFit} function. For additional details on usage, and model
#' assumptions, see the included vignette.
#' 
#'    
#' @param Total the number of individuals in the sample
#' @param J the number of variables observed on each individual
#' @param Rj vector of length J specifying the number of repeated measurements
#'  on each variable
#' @param Nijr an array of dimension (Total, J, max(Rj)) indicating the number
#'  of ranking levels for each replication. For multinomial and bernoulli
#'  variables, Nijr[i,j,r] = 1. For rank variables, Nijr[i,j,r] indicates the
#'  number of candidates ranked.
#' @param K the number of sub-populations
#' @param Vj vector of length J specifying the number of possible candidates
#'  for each variable. For a bernoulli variable Vj[j] = 1. 
#' @param alpha vector of length K which is the hyperparameter for Dirichlet
#'  membership distribution
#' @param theta array of dimension (J,K,max(Vj)) which governs the variable
#'  distributions. theta[j,k,] is the parameter for how sub-population k responds
#'  to the variable j. Initially, any element of theta corresponding to a valid candidate/category must be positive (ie cannot be numerically 0).
#'  However, if the number of candidates differs across variables, any unusued portions of theta should be 0.
#' @param phi array of dimension (Total,K) which specifies the variational
#'  parameters for the membership vectors, lambda. If left NULL, it is initialized 
#'  to a uniformly across all groups. It is highly recommended to leave these NULL.
#' @param delta array of dimension (Total,J,max(Rj), max(N), K) which specifies
#'  the variational parameters for the context vectors Z. If left blank, it is
#'   initialized to a uniformly across all sub-populations. It is highly recommended to leave these NULL.
#' @param dist vector of length J specifying variable types. Options are
#'  "bernoulli", "multinomial" or "rank" corresponing to the distributions
#'   of the observed variables
#' @param obs an array of dimension (Total,J,max(Rj), max(N)) corresponding to 
#' the observed data. For bernoulli random variables, this consists of 0/1's. 
#' For multinomial or rank data this consists of integers 0,1,\ldots,(Vj[j]-1).
#' @return The \code{mixedMemModel} object
#' @examples
#' set.seed(123)
#' Total <- 50 # 50 Individuals
#' J <- 3 # 3 different variables
#' # distributions of each variable
#' dist <- c("multinomial", "bernoulli", "rank") 
#' # 100 repeated measures of the multinomial, 5 repeated measures of the
#' # bernoulli, 1 repeated measure of the rank
#' Rj <- c(100, 5, 1) 
#' 
#' K <- 4 # 4 sub-populations
#' alpha <- rep(.5, K) #hyperparameter for dirichlet distribution
#' 
#' # Number of choices for each variable. Note the Bernoulli must have 1 choice
#' Vj <- c(10, 1, 4) 
#' 
#' 
#' theta <- array(0, dim = c(J, K, max(Vj)))
#' # Parameters governing multinomial
#' theta[1,,] <- gtools::rdirichlet(K, rep(.3, Vj[1]))
#' #parameters governing bernoulli
#' theta[2,,] <- cbind(rbeta(K, 1,1), matrix(0, nrow = K, ncol = Vj[1]-1))
#' theta[3,,] <- cbind(gtools::rdirichlet(K, rep(.3, Vj[3])),
#'  matrix(0, nrow = K, ncol = Vj[1]-Vj[3]))
#' 
#' # generate group memberships
#' lambda <- gtools::rdirichlet(Total, rep(.2,K))
#' 
#' # Candidates selected for each observation. For Multinomial and Bernoulli, this is always 1
#' # For rank data, this will be the number of candidates ranked
#' Nijr = array(0, dim = c(Total, J, max(Rj)))
#' Nijr[,1,c(1:Rj[1])] = 1 # N_ijr is 1 for multinomial variables
#' Nijr[,2,c(1:Rj[2])] = 1 # N_ijr is 1 for Bernoulli variables
#' Nijr[,3,c(1:Rj[3])] = sample.int(Vj[3], size = Total, replace = TRUE)
#' 
#' ## Generate Observations
#' obs = array(0, dim = c(Total, J, max(Rj), max(Nijr)))
#' for(i in 1:Total)
#' {
#'   for(j in 1:J)
#'   {
#'     for(r in 1:Rj[j])
#'     {
#'     # sub-population which governs the observed response
#'       sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) 
#'       if(dist[j] =="multinomial") {
#'        # must be in 0:(Vj[j]-1) so we subtract 1
#'         obs[i,j,r,1] <- sample.int(Vj[j], size = 1, prob = theta[j,sub.pop,c(1:Vj[j])])-1
#'       }
#'       if(dist[j] == "rank") {
#'         # must be in 0:(Vj[j]-1) so we subtract 1
#'         obs[i,j,r,c(1:Nijr[i,j,r])] <- sample.int(Vj[j], size = Nijr[i,j,r],
#'          replace = FALSE, prob = theta[j,sub.pop,c(1:Vj[j])])-1
#'       }
#'       if(dist[j] == "bernoulli"){
#'         obs[i,j,r,1] <- rbinom(n = 1, size = 1, prob = theta[j,sub.pop,1])
#'       }
#'     }
#'   }
#' } 
#' 
#' ## Initialize a mixedMemModel object
#' test_model <- mixedMemModel(Total = Total, J = J,Rj = Rj,
#'  Nijr= Nijr, K = K, Vj = Vj,dist = dist, obs = obs,
#'   alpha = alpha, theta = theta)
#' # Look at Summary of the initialized model
#' summary(test_model)
#' # Plot the current values for theta
#' plot(test_model) 
#' @export

mixedMemModel = function(Total, J, Rj, Nijr, K, Vj, alpha, theta, phi = NULL, delta = NULL, dist, obs)
{
  # Checks if model defaults are used and fills in defaults
  if(is.null(alpha))
  {alpha = rep(1/K,K)}
  
  if(is.null(theta))
  {
    theta = array(0, dim = c(J,K,max(Vj)))
    for(j in 1:J)
    {
      if(dist[j] != "bernoulli")
      {
        for(k in 1:K)
        {
          theta[j,k,] = c(gtools::rdirichlet(1,rep(1,Vj[j])), rep(0, max(Vj)-Vj[j]))
        }
      } else {
        theta[j,,1] =rbeta(K,1,1)
      }
    }
  }
  
  if(is.null(delta))
  {
    delta = array(0, dim = c(Total,J, max(Rj), max(Nijr), K))
    for(i in 1:Total)
    {
      for(j in 1:J)
      {
        for(r in 1:Rj[j])
        {
          for(n in 1:Nijr[i,j,r])
          {
            delta[i,j,r,n,] = rep(1/K,K)
          }
        }
      }
    }
  }
  
  if(is.null(phi))
  {
    phi = array(1/K, dim = c(Total,K))
  }
  #put objects in a list
  model_obj = list(Total, J, Rj, Nijr, K, Vj, alpha, theta, phi, delta, dist, obs);
  names(model_obj) = c("Total", "J", "Rj", "Nijr", "K", "Vj", "alpha","theta", "phi", "delta", "dist" ,"obs")
  class(model_obj) = "mixedMemModel"
  
  dimnames(model_obj$theta) <- list(paste("Var", c(1:J)),
                                    paste("Group", c(1:K)),
                                    paste("Cand", c(0:(max(Vj)-1))))
  
  names(model_obj$alpha) <- paste("Group", c(0:(K-1)))
  
  names(model_obj$Vj) <- paste("Var", c(1:J))
  names(model_obj$Rj) <- paste("Var", c(1:J))
  names(model_obj$dist) <- paste("Var", c(1:J))
  
  dimnames(model_obj$phi) <- list(paste("Ind", c(1:Total)),
                                  paste("Group", c(0:(K-1))))
  
  dimnames(model_obj$delta) <- list(paste("Ind", c(1:Total)),
                                    paste("Var", c(1:J)),
                                    paste("Rep", c(1:max(Rj))),
                                    paste("Rank", c(1:max(Nijr))),
                                    paste("Group", c(0:(K-1))))
  
  dimnames(model_obj$Nijr) <- list(paste("Ind", c(1:Total)),
                                   paste("Var", c(1:J)),
                                   paste("Rep", c(1:max(Rj))))
  
  dimnames(model_obj$obs) <- list(paste("Ind", c(1:Total)),
                                  paste("Var", c(1:J)),
                                  paste("Rep", c(1:max(Rj))),
                                  paste("Rank", c(1:max(Nijr))))
  #check for valid model parameters
  checkModel(model_obj)
  return(model_obj)
}


#' Summary of a Mixed Membership Model
#' 
#' Summary of a mixedMemModel object 
#' 
#' Prints summary data for a mixedMemModel object which includes the
#' ELBO, the dimensions of the model and details about each variable.
#'  
#' @param object the mixedMemModel object to be summarized
#' @param ... additional parameters
#' @seealso mixedMemModel
#' @export
summary.mixedMemModel = function(object,...)
{
  cat("== Summary for Mixed Membership Model ==\n")
  cat(paste("Total: ", object$Total, "\t\t K: ",object$K, "\t\t ELBO: ",round(computeELBO(object),2),"\n\n" ,sep = ""))
  df = data.frame(paste("  ",c(1:object$J),sep = ""), object$dist, object$Rj,
                  object$Vj)
  colnames(df) = c("Variable  ", "Variable Type    ", "Replicates    ", "Categories  ")
  print(df, row.names = FALSE, right = FALSE)
}

#' Plot a Mixed Membership Model
#' 
#' Visual representation of a mixedMemModel object 
#' 
#' @param x the mixedMemModel object to be plotted
#' @param type which estimated parameters to plot; either "theta" or "membership" 
#' @param compare model to compare. For type = "theta", compare should be an array the same size as x$theta
#' for type = "membership", compare should be a matrix the same size as x$phi
#' @param main title for chart
#' @param varNames vector of names for each variable if plotting theta
#' @param groupNames vector of labels for each sub-population
#' @param nrow number of rows for the grid of plots
#' @param ncol number of columns for the grid of plots. if plotting theta, this
#' must be K, if plotting membership, this can be specified
#' @param indices when plotting memberships, which individuals to plot. When plotting theta, which variables to plot 
#' @param fitNames vector of labels for each fit
#' @param ... additional parameters
#' @seealso mixedMemModel, vizTheta, vizMem
#' @export
plot.mixedMemModel = function(x, type = "theta" , compare = NULL,
                              main = NULL,
                              varNames = NULL,
                              groupNames = NULL,
                              nrow = NULL, ncol = NULL, indices = NULL, fitNames = NULL,...){
  if(type =="theta") {
    if(is.null(main)){
      main = "Estimated Theta"
    }
    vizTheta(x, compare = compare, main = main, varNames = varNames,
             groupNames = groupNames, nrow = nrow, fitNames = fitNames, indices = indices)
  } else if (type == "membership") {
    if(is.null(main)){
      main = "Estimated Memberships"
    }
    vizMem(x,compare = compare, main = main, nrow = nrow, ncol = ncol,
           indices = indices, fitNames = fitNames, groupNames = groupNames)
  }
}
