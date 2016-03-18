#' Mixed Membership Post-Processing
#'    
#' \code{findLabels} finds the optimal permutation of labels that minimizes 
#' the weighted squared difference between the arrays of subpopulation parameters from a fitted mixed membership
#' model, \eqn{\theta} and a given comparison model. 
#' 
#' 
#' Mixed Membership models are invariant to permutations of the sub-population labels; swapping the names of each sub-population yields an equivalent model. 
#' The ordering of the labels in a fitted model is often dependent on the initialization points.
#' The function \code{findLabels} selects a permutation of the sub-population labels that best matches a given comparison model 
#' by minimizing squared differences between the \eqn{\theta} arrays. 
#' 
#' \eqn{Loss = \sum_j \sum_k [\sum_v (\hat\theta_{k,v} - \theta_{k,v})^2]}
#' 
#' If K, number of sub-populations, is small, the method searches through all K! permutations of the sub-population labels and 
#' select the permutation which minimizes the loss. If K is large, a greedy algorithim can be used instead. This
#' algorithm selects the best match for each fitted sub-population starting with the group with the largest fitted 
#' relative frequency.
#'  
#' @param model the fitted \code{mixedMemModelMCMC} or \code{mixedMemModelVI} object.
#' @param comparison an array of the same dimensions as model$theta which contains the subpopulation parameters from another model.
#'  \code{findLabels} will return a permutation of the labels of \code{model} which match to \code{comparison} most closely.
#' @param exhaustive a boolean for whether an exhaustive search should be performed. If false, a greedy algorithim is used instead.
#' @return \code{findLabels} returns a list with two objects: \code{perm} and \code{loss}. \code{perm} is the optimal permutation of the labels with respect to the squared error loss.
#' \code{loss} is the calculated value of the weighted squared error loss (shown above) for the optimal permutation.
#' @seealso permuteLabels
#' @examples
#' 
#' \dontrun{
#' # See mixedMemModelMCMC or mixedMemModelVI documentation for how to generate data and instantiate objects
#' # After the data as been generated, we initialize the array of sub-population parameters (theta) 
#' # according to a permutation of the true labeling
#' set.seed(123)
#' perm <- sample.int(K, size = K, replace = FALSE)
#' theta.perm <- theta_truth[,perm,]
#' test_model <- mixedMemModelVI(Total = Total, J = J,Rj = Rj, Nijr= Nijr,
#'  K = K, Vj = Vj,dist = dist, obs = obs, alpha = alpha, theta = theta.perm)
#' out <- mmVIFit(test_model)
#' opt.perm <- findLabels(out, theta_truth)
#' opt.perm
#' 
#' # produce mixedMemModel object with sub-population labels permuted to best match
#' # the comparison model
#' out <- permuteLabels(out, opt.perm$perm)
#' }
#' @export
findLabels <- function(model, comparison,  exhaustive = FALSE) {
  
    fitted.set <- model$theta
    K <- model$K
 
    optimal.perm = c(1:K)
    
    
    loss <- 0
    if(exhaustive) {
      # list all permutations of K
      perms <- gtools::permutations(K, K)
      # get loss of 1st permutation
      loss <- sum((fitted.set[, perms[1, ], , drop = F] - comparison)^2)
      
      # check each possible permutation 
      for(i in 2:factorial(K)) {
        # calculate squared error loss for permutation i
        loss.i <- sum((fitted.set[,perms[i, ], , drop = F] - comparison)^2)
        
        # check if squared error loss is less than current best  
        if(loss.i < loss) {
            loss <- loss.i
            optimal.perm <- perms[i,] 
        }

      }
      
    
    } else {   # Greedy search 
      
      # get search ordering
      if(class(model) == "mixedMemModelMCMC"){
        search.order <- order(model$ksi, decreasing = T)
      } else if (class(model) == "mixedMemModelVI") {
        search.order <- order(model$alpha, decreasing = T)
      }

      selected <- c()
      loss.j <- rep(0,K)
      
      for(i in 1:K) {
        for(j in 1:K) {
          loss.j[j] <- sum((fitted.set[, search.order[i] , , drop = F] - comparison[,j, , drop = F])^2) 
        }
        
        # remove already selected groups
        loss.j[selected] <- max(loss.j) + 1
        selected <- c(selected, which.min(loss.j))
        
        optimal.perm[ search.order[i] ] <- which.min(loss.j)
        loss <- loss + min(loss.j)
      }
    }
    return(list(perm = optimal.perm, loss = loss))
}

#' Mixed Membership Post-Processing
#' 
#' Mixed Membership models are invariant to permutations of the sub-population labels; swapping the names of each sub-population yields an equivalent model. 
#' The ordering of the labels in a fitted model is dependent on the initialization points of the variational EM algorithim.
#' The \code{permuteLabels} function returns a \code{mixedMemModel} object where the labels (for \eqn{\theta}, \eqn{\phi}, \eqn{\delta} and \eqn{\alpha}) have been permuted
#' according a given permutation of the integers 1 through K. The \code{findLabels} function can be used to find a permutation of the labels which
#' most closely matches another fitted model. 
#' 
#' @param model a fitted \code{mixedMemModel} object which will be relabeled.
#' @param perm a vector of length K with integers 1:K. This is the permutation by which to relabel the \code{mixedMemModel} object such that
#' group i in the returned mixedMemModel object corresponds to group \code{perm}[i] from the input mixedMemModel object.
#' @return \code{permuteLabels} returns a \code{mixedMemModel} object such that
#' group i in the returned \code{mixedMemModel} object corresponds to group perm[i] from the input \code{mixedMemModel} object
#' @seealso findLabels
#' @export
permuteLabels <- function(model, perm) {
  
  if(length(perm)!= model$K) {
    stop("Error: perm  must be of length model$K")
  }
  
  
  if(class(model) == "mixedMemModelVI"){
    # Permute relevant quantities for mixedMemModelVI
    theta <- model$theta[, perm, , drop = F]
    alpha <- model$alpha[perm]
    
    phi <- model$phi[, perm, drop = F]
    delta <- model$delta[, , , , perm, drop = F]
    
    out <- mixedMemModelVI(Total = model$Total, J = model$J, Rj = model$Rj, Nijr = model$Nijr,
                           K = model$K, Vj = model$Vj, alpha = alpha, theta = theta,
                           phi = phi, delta = delta, dist = model$dist, obs = model$obs)

  } else if (class(model) == "mixedMemModelMCMC"){
    # Permute relevant quantities for mixedMemModelMCMC
    
    theta <- model$theta[, perm, , drop = F]
    ksi <- model$ksi[perm]
    lambda <- model$lambda[, perm]
    tau <- model$tau[, perm, , drop = F]
    
    Z <- array(0, dim = dim(model$Z))
    
    for(k in 1:model$K){
      Z <- Z + ifelse(Z == k, perm[k], 0)
    }
    
     out <- mixedMemModelMCMC(Total = model$Total, J = model$J, Rj = model$Rj, 
                              K = model$K, Vj = model$Vj, dist = model$dist, obs = model$obs,
                              theta = theta, lambda = lambda, Z = Z, alpha = model$alpha,
                              ksi = ksi, phi = model$phi, delta = model$delta, 
                              phi = model$phi, tau = tau, P = model$P, fixedObs = model$fixedObs, extended = !is.null(model$P))
    
  }

  return(out)
}