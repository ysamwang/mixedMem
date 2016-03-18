#' Mixed Membership Visualization
#' 
#' 
#' \code{entropyPlot} plots a histogram of the entropy of each individuals estimated mixed membership. The entropy is
#' is roughly the number of profiles needed to model each individual (Gormley and Murphy 2014; White et al 2012).
#' 
#' @param model the \code{mixedMemModel} the model to be plotted
#' @param col the color of the bars
#' @param main the main title
#' @param xlab the label for the x-axis
#' @export
entropyPlot = function(model, col = "white", main = "Exponentiated Entropy",
                       xlab = "Effective Profiles", ...) {
  if(class(model) == "mixeMemModelVI"){
    mem <- model$phi / rowSums(model$phi)
  } else if (class(model) == "mixedMemModelMCMC"){
    mem <- model$lambda
  }
  
  entropy <- apply(mem, MAR = 1, function(x){return( exp(-sum(log(x) * x) )) } )
  
  hist(entropy, main = main, col = col, xlab = xlab, ...)
}
