#' Fit Mixed Membership models using a Metropolis Hastings within Gibbs MCMC sampler
#' 
#' Takes samples form the posterior distribution of a mixed membership model using a 
#' Metropolis-Hastings within Gibbs sampler. 
#' 
#' \code{mmMCMCFit} draws samples from the posterior distribution of a mixed membership model
#' given a set of observed data. The sampler takes Metropolis Hastings steps to sample the \eqn{\alpha_0}
#' and \eqn{\xi} parameters and samples the remaining \eqn{\theta}, \eqn{\lambda_i} and \eqn{Z} parameters
#' with a Gibbs sampler.   
#' 
#' @param model a \code{mixedMemModelMCMC} object created by the \code{mixedMemModelMCMC} constructor
#' @param burnIn non-negative integer indicating the number of burn in samples before recording the first sample
#' @param samples positive integer indicating the number of samples to record
#' @param thin positive integer indicating how to thin the samples
#' @param print positive integer indicating how often to print an update to the R console
#' @param fileNames list of files names to write samples
#' @param newFiles 0 if samples should be appended to existing files; 1 if samples should overwrite any existing files 
#' @param omega tuning parameter for MH step for alpha.
#' @param eta tuning parameter for MH step for ksi
#' @param whichWrite which parameters to write to disk. Writing the individual parameters can significantly increase computational
#' time when the number of individuals is large
#' @param extended boolean whether to estimate the extended model or not
#' @seealso mixedMemModelMCMC
#' @export
mmMCMCFit <- function(model, burnIn = 20000, samples = 1000, thin = 10, print = 100, fileNames, newFiles = 1,
                      omega = 100, eta = 1, whichWrite = c(1, 1, 1, 0, 0, 1, 0)) {
  
  # check to make sure the model is of the right class
  if(class(model) != "mixedMemModelMCMC"){
    stop("Input model must be of class mixedMemModelMCMC")
  }
  
  # checks inputs for internal consistency
  checkModelMCMC(model) 
  
  # In C code, Z must be 0:(K-1), but in R it is it 1:K so we modify it before inputting to the c++ function
  model$Z <- model$Z -1
  mcmcInputC(model = model, burnIn = burnIn, samples = samples, thin = thin, print = print, fileNames = fileNames,
                      newFiles = newFiles, omega = omega, eta = eta, whichWrite = whichWrite)
  
  
}
