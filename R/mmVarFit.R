#' Fit  Mixed Membership models using variational EM
#' 
#' Main function of the \code{mixedMem} package. Fits parameters \eqn{\phi} and \eqn{\delta} for the variational
#'  distribution of latent variables as well as psuedo-MLE estimates for 
#'  the population hyperparameters \eqn{\alpha} and \eqn{\theta}. See documentation for
#'  \code{mixedMemModel} or the package vignette for a more detailed description of 
#'  the variables/notation in a mixed membership model. 
#' 
#'  
#' \code{mmVarFit} selects psuedo-MLE estimates for \eqn{alpha} and \eqn{theta} and approximates
#' the posterior distribution for the latent variables through a mean field variational approach.
#' Using Jensen's inequality, we derive the lower bound on the RHS (sometimes called the ELBO) for the log-likelihood
#' of the data.
#' 
#' \eqn{P(obs |\alpha, \theta) \ge E_Q{\log[p(X,Z, \Lambda)]} - E_Q{\log[Q(Z, \Lambda|\phi, \delta)]}}
#' 
#' where
#' 
#' \eqn{Q = \prod_i [Dirichlet(\lambda_i|\phi_i) \prod_j^J \prod_r^{R_j} \prod_n^{N_{i,j,r}}Multinomial(Z_{i,j,r,n}|\delta_{i,j,r,n})]}
#' 
#' It can be shown that maximizing the ELBO with respect to \eqn{\phi} and \eqn{\delta}   
#' minimizes the KL divergence, between a tractable variational distribution and the true posterior.
#' We can simultaneously pick psuedo-MLE hyperparameters \eqn{\alpha} and \eqn{\theta} to maximize the lower bound on 
#' the log-likelihood of the observed data.
#'
#' The method uses an EM approach. The E step considers the hyperparameters fixed and picks appropriate variational
#' parameters to minimize the KL divergence. On the M step
#' the variational parameters are fixed and hyperparmaters are selected which maximize the lower bound 
#' on the log-likelihood. Aside from printStatus, printMod, and stepType, the authors do not recommend changing the other
#' parameters.
#'   
#' @references
#' Beal, Matthew James. Variational algorithms for approximate Bayesian inference. Diss. University of London, 2003.
#'
#' @param model a \code{mixedMemModel} object created by the \code{mixedMemModel} constructor
#' @param printStatus Options are 1 or 0. 1 will print status updates, 0 will not print status updates.
#' @param printMod Positive integer which specifies how often to print status updates. The status will be printed at each step which is a multiple of printMod.
#' @param stepType An integer (0,1,2,3)  specifying what steps to fit. 0 only performs an E-Step; this can be used
#' to find the held out ELBO. 1 performs E-steps and fits theta but keeps alpha constant. 2 performs E-steps and fits alpha, but keeps theta constant. 3 completes full
#' E and M steps.
#' @param maxTotalIter The maximum total steps before termination. A full E and M step together count as 1 step. If this maximum is ever achieved, a warning message will be printed at convergence.
#' @param maxEIter The maximum iterations for each the E-Step. If this maximum is ever achieved, a warning message will be printed at convergence.
#' @param maxAlphaIter The maximum iterations when fitting alpha. If this maximum is ever achieved, a warning message will be printed at convergence.
#' @param maxThetaIter The maximum iterations when fitting theta. If this maximum is ever achieved, a warning message will be printed at convergence.
#' @param maxLSIter The maximum backtracking iterations in the line search for updating alpha and theta for rank data 
#' @param elboTol The convergence criteria for the EM Algorithim. When the relative increase in the ELBO is less than the convergence critiera,
#' the algorithim converges 
#' @param alphaTol The convergence criteria for updates to alpha. When the relative increase in the ELBO is less than the convergence critiera,
#' the update to alpha converges
#' @param thetaTol The convergence criteria for updates to theta. When the relative increase in the ELBO is less than the convergence critiera,
#' the update to theta converges
#' @param aNaught The first step size in the backtracking line search used to update alpha and theta for rank data
#' @param tau The backtracking factor in the backtracking line search used to update alpha and theta for rank data
#' @param bMax The number of iterations for the interior point method for fitting theta for rank data
#' @param bNaught The initial scaling factor in the interior point method for fitting theta for rank data
#' @param bMult The factor by which bNaught is multiplied by in each iteration of the interior point methodfor fitting theta for rank data
#' @param vCutoff The cutoff for Vj at which a gradient ascent method is used instead of the Newton Raphson interior point method. This is used to avoid inverting
#' a large matrix 
#' @param holdConst an vector of integers specifying groups to be held constant during the estimation procedure. The estimation algorithim will hold the
#'  theta parameters of these specific groups constant, but update all other parameters. The group numbers range from 0 to K - 1.
#' @return a \code{mixedMemModel} containing updated variational parameters and hyperparameters
#' @seealso mixedMemModel
#' @examples
#' 
#' ## Generate Data 
#' Total <- 30 #30 Individuals
#' J <- 2 # 2 variables
#' dist <- rep("multinomial",J) #both variables are multinomial
#' Rj <- rep(100,J) #100 repititions for each variable
#' #Nijr will always be 1 for multinomials and bernoulli's
#' Nijr <- array(1, dim = c(Total, J, max(Rj))) 
#' K <- 4 # 4 sub-populations
#' alpha <- rep(.5, K) #hyperparameter for dirichlet distribution
#' Vj <- rep(5, J) #each multinomial has 5 options
#' theta <- array(0, dim = c(J, K, max(Vj)))
#' theta[1,,] <- gtools::rdirichlet(K, rep(.3, 5))
#' theta[2,,] <- gtools::rdirichlet(K, rep(.3, 5))
#' lambda <- gtools::rdirichlet(Total, rep(.6,K))
#' obs = array(0, dim = c(Total, J, max(Rj), max(Nijr)))
#' for(i in 1:Total)
#' {
#'  for(j in 1:J)
#'  {
#'    for(r in 1:Rj[j])
#'    {
#'      for(n in Nijr[i,j,r])
#'      {
#'      # sub-population which governs the multinomial
#'      sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) 
#'      #Note that observations must be from 0:(Vj-1)
#'      obs[i,j,r,n] <- sample.int(Vj[j], size = 1,prob = theta[j,sub.pop,])-1       }
#'    }
#'  } 
#' }
#' 
#' ## Initialize a mixedMemModel object
#' test_model <- mixedMemModel(Total = Total, J = J,Rj = Rj, Nijr= Nijr,
#'  K = K, Vj = Vj,dist = dist, obs = obs,
#'   alpha = alpha, theta = theta+0)
#' 
#' ## Fit the mixed membership model
#' out <-mmVarFit(test_model)
#' @export
mmVarFit = function(model, printStatus = 1,
                    printMod = 1, stepType = 3,
                    maxTotalIter = 500, maxEIter = 1000,
                    maxAlphaIter = 200, maxThetaIter = 1000,
                    maxLSIter = 400, elboTol = 1e-6, alphaTol = 1e-6,
                    thetaTol = 1e-10, aNaught = 1.0, tau = .899,
                    bMax = 3, bNaught = 1000.0, bMult = 1000.0, vCutoff = 13, holdConst = c(-1)) {
  output = model
  names(output) = c("Total", "J", "Rj", "Nijr", "K", "Vj", "alpha","theta", "phi", "delta", "dist" ,"obs")
  output$alpha = (model$alpha+0)
  output$theta = (model$theta+0)
  output$phi = (model$phi+0)
  output$delta = (model$delta+0)
  
  checkModel(output) # R function which checks inputs
  print("Model Check: Ok!")
  print("<== Beginning Model Fit! ==>")

  varInfInputC(output, printStatus, printMod, stepType, maxTotalIter, maxEIter, maxAlphaIter,
               maxThetaIter, maxLSIter, elboTol, alphaTol, thetaTol, aNaught, tau, bMax, bNaught, 
               bMult, vCutoff, holdConst) # R wrapper function
  return(output)
}



#' Compute a lower bound on the log-likelihood (ELBO)
#' 
#' Computes the value of a lower bound on the log-likelihood, also called the ELBO, for a mixed membership model.
#'
#' The lower bound (ELBO) is the objective function in the variational EM algorithim. It is a function of the latent variables (\eqn{\phi} and \eqn{\delta}) and the hyperparameters (\eqn{\alpha} and \eqn{\theta}) and
#' can be derived by Jensen's inequality.
#'      
#' \eqn{P(obs |\alpha, \theta) \ge E_Q{\log[p(X,Z, \Lambda)]} - E_Q{\log[Q(Z, \Lambda|\phi, \delta)]}}
#'    
#' @param model a \code{mixedMemModel} object created by the \code{mixedMemModel} constructor
#' @return value of the lower bound on the log-likelihood
#' @export
computeELBO = function(model)
{
  checkModel(model)
  return(computeElboC(model))
}