#' Mixed Membership Post-Processing
#'    
#' Finds the permutation of labels that minimizes 
#' the weighted squared error loss between the fitted theta and a comparison model. 
#' 
#' 
#' Mixed Membership models are invariant to permutations of the sub-population labels and the ordering of the labels in a fitted model
#' is dependent on the initialization points of the variational EM algorithim. The \code{findLabels} function selects an 
#' optimal permutation of the labels to match a given comparison model.
#' The loss function is the weighted sum of squared differences where the weights are determined by the relative frequency of each group.  
#' 
#' \eqn{Loss = \sum_j \sum_k \alpha_k/\alpha_0 [\sum_v (\hat\theta_{k,v} - \theta_{k,v})^2]}
#' where \eqn{\alpha_0 = \sum_k \alpha_k}
#' 
#' If K, number of sub-populations, is small, the method can search through all K! permutations and 
#' select the permutation which minimizes the loss. If K is large, a greedy algorithim can be used instead. This
#' algorithim selects the best match for each fitted sub-population starting with the group with the largest fitted 
#' relative frequency.
#'  
#' @param model the fitted \code{mixedMemModel} object
#' @param training the comparison theta. This should be an array the same dimensions as model$theta
#' @param exhaustive a boolean for whether an exhaustive search should be performed. If false, a greedy algorithim is used instead
#' @return perm optimal permutation of the labels with respect to squared error loss
#' @return loss the sum of squared error loss of the optimal permutation weighted by relative frequency
#' @seealso permuteLabels
#' @examples
#' 
#' \dontrun{
#' # See mixedMemModel documentation for how to generate data and instantiate a mixedMemModel object
#' # After the data as been generated, we initialize theta to a permutation of the true labeling
#' set.seed(123)
#' perm = sample.int(K, size = K, replace = FALSE)
#' theta1 = theta[,perm,]
#' test_model <- mixedMemModel(Total = Total, J = J,Rj = Rj, Nijr= Nijr, K = K, Vj = Vj,dist = dist,
#'  obs = obs, alpha = alpha, theta = theta1)
#' out <- mmVarFit(test_model)
#' opt.perm <- findLabels(out, theta)
#' opt.perm
#' 
#' # produce mixedMemModel object with labels permuted to match ground truth
#' out = permuteLabels(out, opt.perm$perm)
#' }
#' @export
findLabels = function(model, training,  exhaustive = TRUE)
{
    fitted.set = model$theta
    K = model$K
 
    optimal.perm = c(1:K)
    weight = model$alpha/sum(model$alpha)
    loss = 0
    if(exhaustive)
    {
    perms = gtools::permutations(K,K)
    loss = sum((fitted.set[,perms[1,],]-training[,c(1:K),])^2)
    
    for(i in 2:factorial(K))
      {
        if(dim(training)[1]==1) {
          diff = (fitted.set[,perms[i,],]-training[,c(1:K),])^2
        } else {
          diff = aperm((fitted.set[,perms[i,],]-training)^2, c(2,1,3))
        }
        loss.i = sum(diff*weight)
        if(loss.i<loss)
        {
          loss = loss.i
          optimal.perm = perms[i,] 
        }
      }
    } else {
      priority = order(weight, decreasing = T)
      selected = c()
      loss.j = rep(0,K)
      for(i in 1:K)
      {
        for(j in 1:K)
        {
          loss.j[j] = sum((fitted.set[,priority[i],] - training[,j,])^2)*weight[i] 
        }
        loss.j[selected] = max(loss.j)+1
        selected = c(selected, which.min(loss.j))
        optimal.perm[priority[i]] = which.min(loss.j)
        loss = loss + min(loss.j)
      }
    }
    return(list(perm = optimal.perm, loss = loss))
}

#' Mixed Membership Post-Processing
#' 
#' Mixed Membership models are invariant to permutations of the sub-population labels and the ordering of the labels
#' is dependent on the initialization of the variational EM algorithim.
#' Given a permutation, which can be found using \code{findLabels},
#' the method returns a \code{mixedMemModel} object with permuted labels
#' 
#' @param model the fitted \code{mixedMemModel} object
#' @param perm the permutation by which to relabel the \code{mixedMemModel} object. Must be a vector of length model$K.
#' @return a permuted \code{mixedMemModel} object
#' @seealso findLabels
#' @export

permuteLabels = function(model, perm)
{
  if(length(perm)!= model$K)
  {stop("Error: perm  must be of length model$K")}
  
  out = model
  
  out$alpha = out$alpha[perm]
  names(out$alpha) <- names(model$alpha)
  for(j in 1:out$J)
  {
    out$theta[j,,] = out$theta[j,perm,]
  }
  dimnames(out$theta) <- dimnames(model$theta)
  
  out$phi = out$phi[,perm]
  for(i in 1:out$Total)
  {
  for(j in 1:out$J)
  {
    for(r in 1:out$Rj[j])
    {
      for(n in 1:out$Nijr[i,j,r])
      {
        out$delta[i,j,r,n,] = out$delta[i,j,r,n,perm]
      }
    }
  }
  }
  
  dimnames(out$delta) <- dimnames(model$delta)
  dimnames(out$phi) <- dimnames(model$phi)
  
  return(out)
}