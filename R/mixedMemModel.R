#' Constructor for a Mixed Membership Model Object
#' 
#' Constructor for a \code{mixedMemModel} object which can be used for analysis 
#' in the mixedMem package.
#' 
#' The function returns an object of \code{mixedMemModel} class. This object contains dimensions of the model,
#' the observed data, and the estimates of the model parameters. Once a \code{mixedMemModel} object is created,
#' the specified model can be fit for the data using the \code{mmVarFit} function. For additional details on usage, and model
#' assumptions, see the corresponding vignette "Fitting Mixed Membership Models using \texttt{mixedMem}".
#' 
#'    
#' @param Total the number of individuals in the sample
#' @param J the number of variables observed on each individual
#' @param Rj vector of length J specifying the number of repeated measurements
#'  for each variable
#' @param Nijr an array of dimension (Total, J, max(Rj)) indicating the number
#'  of ranking levels for each variable and each replication. For multinomial and bernoulli
#'  variables, Nijr[i,j,r] = 1. For rank variables, Nijr[i,j,r] indicates the
#'  number of alternatives ranked.
#' @param K the number of sub-populations
#' @param Vj vector of length J specifying the number of possible candidates
#'  for each variable. For a bernoulli variable Vj[j] = 1. 
#' @param alpha vector of length K which is the hyperparameter for Dirichlet
#'  membership distribution
#' @param theta array of dimension (J,K,max(Vj)) which governs the variable
#'  distributions. theta[j,k,] is the parameter for how sub-population k responds
#'  to the variable j. For instance, if variable j is a Bernoulli variable, theta[j,k,1] is the probability of success; if
#'  variable j is a multinomial variable, theta[j,k, 1:Vj[j]] is the probability for each of the Vj[j] categories ; if variable j
#'  is a rank variable, theta[j,k, 1:Vj[j]] are the support parameters for each of the Vj[j] alternatives. Since the dimension of the relevant parameters
#'  can differ across variables, any unused elements of the array should be set to 0, while all other elements should be positive.
#' @param phi array of dimension (Total,K) which specifies the variational
#'  parameters for the membership vectors, lambda. The default group membership initialization is uniform across all groups (phi[i,k] = 1/K for all k).
#'  The default initialization is highly recommended.
#' @param delta array of dimension (Total,J,max(Rj), max(N), K) which specifies
#'  the variational parameters for the context vectors Z. The default initialization is
#'   uniform across all sub-populations (delta[i, j, r, n, k] = 1/K for all k). The default initialization is highly recommended.
#' @param dist vector of length J specifying variable types. Options are
#'  "bernoulli", "multinomial" or "rank" corresponing to the distributions
#'   of the observed variables
#' @param obs an array of dimension (Total,J,max(Rj), max(N)) corresponding to 
#' the observed data. For bernoulli random variables, the data consist of 0/1's. 
#' For multinomial or rank data the data consist of integers 0,1,\ldots,(Vj[j]-1).
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
#' # Number of categories/alternatives for each variable. Note the Bernoulli must Vj = 1
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
#' # Alternatives selected for each observation. For Multinomial and Bernoulli, this is always 1
#' # For rank data, this will be the number of candidates ranked
#' Nijr = array(0, dim = c(Total, J, max(Rj)))
#' Nijr[,1,c(1:Rj[1])] = 1 # N_ijr is 1 for multinomial variables
#' Nijr[,2,c(1:Rj[2])] = 1 # N_ijr is 1 for Bernoulli variables
#' Nijr[,3,c(1:Rj[3])] = sample.int(Vj[3], size = Total, replace = TRUE)
#' 
#' # generate random sample of observations
#' sampleMixedMem <- rmixedMem(Total, J, Rj, Nijr, K, Vj,
#' dist, theta, alpha)
#' 
#' ## Initialize a mixedMemModel object
#' test_model <- mixedMemModel(Total = Total, J = J,Rj = Rj,
#'  Nijr= Nijr, K = K, Vj = Vj,dist = dist, obs = sampleMixedMem$obs,
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
#' Provides a summary of a mixedMemModel object 
#' 
#' Generic S3 function to provide a summary for a given \code{mixedMemModel} object. The function 
#' prints the the ELBO, the dimensions of the model and each variable type.
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
#' Generic function to produce visual representation of a mixedMemModel object. This function calls either the \code{vizTheta} or
#' the \code{vizMem} function. 
#' 
#' @param x the mixedMemModel object to be plotted
#' @param type which estimated parameters to plot; either "theta" or "membership". \code{vizTheta} is called when the type is "theta"
#' and \code{vizMem} is called when the type is "membership". 
#' @param compare object for comparison. When type = "theta", \code{compare} should be an array the same size as x$theta.
#' When type = "membership", \code{compare} should be a matrix the same size as x$phi.
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
