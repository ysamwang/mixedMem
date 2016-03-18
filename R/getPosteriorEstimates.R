#' Mixed Membership Post-Processing for MCMC method
#'    
#' \code{getPosteriorEstimates} 
#' 
#'  
#' @param model the fitted \code{mixedMemModelMCMC} object.
#' @param fileNames a vector of strings containing the file paths for the MCMC output
#' The string should be of length 7 for theta, alpha, ksi, lambda, Z, P and StayerStatus.
#' If the particular parameter was not recorded, insert an empty string ("") as a placeholder 
#' @param whichWrite, a vector of length 7 indicating which files to be read. 1 indicates read,
#' 0 indicates skip
#' @export
getPosteriorEstimates <- function(model, fileNames, whichWrite = c(1, 1, 1, 0, 0, 1, 0)) {
  
  if(model$extended){
    names.list <- c("Theta", "Alpha", "Ksi", "Lambda", "Z", "P", "StayerStatus")
  } else {
    names.list <- c("Theta", "Alpha", "Ksi", "Lambda", "Z")
  }
  
  for(n in 1:length(names.list)){
    if(whichWrite[n] == 1){
      name <- paste("out_", names.list[n], sep = "")
      assign(name, read.csv(fileNames[n], header = F))
    }
  }
  
  
  theta.hat <- NULL
  theta.var <- NULL
  if(whichWrite[1]){
    theta.hat <- array(colMeans(out_Theta), dim = dim(model$theta))
    theta.var <- array(apply(out_Theta, MAR = 2, FUN = var), dim = dim(model$theta))
  }
    
  
  alpha.hat <- alpha.var <- NULL
  if(whichWrite[2]){
    alpha.hat <- mean(unlist(out_Alpha))
    alpha.var <- var(unlist(out_Alpha))
  }
  
  ksi.hat <- NULL
  ksi.var <- NULL
  if(whichWrite[3]){
    ksi.hat <- colMeans(out_Ksi)
    ksi.var <- apply(out_Ksi, MAR = 2, FUN = var)
  }
  
  lambda.hat <- NULL
  lambda.var <- NULL
  if(whichWrite[4]){
    lambda.hat <- matrix(colMeans(out_Lambda), nrow = nrow(model$lambda), nrow = ncol(model$lambda) )
    lambda.var <- matrix(apply(out_Lambda, MAR = 2, FUN = var), nrow = nrow(model$lambda), nrow = ncol(model$lambda) )
  }
  
  # TODO: Think about what summary of Z makes sense
  
  p.hat <- p.var <- NULL
  s.hat <-  NULL
  
  if(model$extended){
    if(whichWrite[6]){
      if(!is.null(model$P)) {
        p.hat <- colMeans(out_P)
        p.var <- apply(out_P, MAR = 2, FUN = var)
      }
    }
    
    if(whichWrite[7]){
        S.hat <- colMeans(out_StayerStatus)
    }
  }
  
  ret <- list(theta.hat = theta.hat, theta.var = theta.var,
              alpha.hat = alpha.hat, alpha.var = alpha.var,
              ksi.hat = ksi.hat, ksi.var = ksi.var,
              lambda.hat = lambda.hat, lambda.var = lambda.var,
              p.hat = p.hat, p.var = p.var,
              s.hat = s.hat)
  # class(ret) <- "mixedMemMCMCOut"
  return(ret)
}




#' Mixed Membership Post-Processing for MCMC method
#'    
#' \code{plotPosterior} 
#' 
#'  
#' @param model the fitted \code{mixedMemModelMCMC} object.
#' @param fileName a string containing the file paths for the MCMC output
#' @param parameter which parameter to plot
#' @param j the jth element to plot
#' @param k the kth element to plot
#' @param v the vth element to plot
#' @param s the sth element to plot
#' @param type histogram or trace plot
#' @param addToExisting whether to plot on existing plot or new plot
#' @param ... additional commands passed to plotting function
#' @export
plotPosterior <- function(model, fileName, parameter,
                                   addToExisting = F, plotType = "trace", j = 1, k = 1, v = 1, s = 1, ...) {
  # whether to add to existing plot
  if(addToExisting == F){
    plottingFunc <- plot
  } else {
    plottingFunc <- lines
  }
  
  parameter <- tolower(parameter)
  
  out <- read.csv(fileName, header = F)
  
  if (plotType == "trace"){
    if(parameter == "alpha"){
      plottingFunc(unlist(out), ...)
    } else if (parameter == "ksi"){
        plottingFunc(out[,k], ...)
    } else if(parameter == "theta") {
      ind <- j + (model$J) * (k-1) + (model$J) * (model$K) * (v - 1) 
      plottingFunc(out[, ind], ...)
    } else if(parameter == "p"){
      plottingFunc(out[, s], ...)
    }
  } else if (plotType == "hist"){
    if(parameter == "alpha"){
      hist(unlist(out), add = addToExisting, ...)
    } else if (parameter == "ksi"){
      hist(out[,k], add = addToExisting, ...)
    } else if(parameter == "theta") {
      ind <- j + (model$J - 1) * k + (model$J - 1) * (model$K - 1) * v 
      hist(out[, ind],  add = addToExisting, ...)
    } else if(parameter == "p"){
      hist(out[, s],  add = addToExisting, ...)
    }
  }
}





