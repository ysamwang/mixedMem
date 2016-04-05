#' Mixed membership models post-processing
#' 
#' 
#' \code{hellingerDistances} computes hellinger distances between each sub-population for each variable
#' 
#' @param model the \code{mixedMemModelVI} object that will be plotted.
hellingerDistances <- function(model){
  distances <- array(0, dim = c(model$J, model$K, model$K))
    for (j in 1:model$J) {
      for (k1 in 1:model$K){
        for (k2 in 1:k1) {
          distances[j, k1, k2] <- sqrt(sum((sqrt(model$theta[j, k1, 1:model$Vj[j]]) -
                                         sqrt(model$theta[j, k2, 1:model$Vj[j]]))^2)) / sqrt(2)
          distances[j, k2, k1] <- distances[j, k1, k2] 
        }
      }
    }

  return(distances)
}