#' Mixed Membership Post-Processing
#'    
#' Finds the permutation of labels that minimizes 
#' the weighted squared error loss between the fitted theta and a ground truth theta. 
#' 
#' 
#' Mixed Membership models are invariant to permutations of the sub-population labels and the ordering of the labels in a fitted model
#' is dependent on the initialization points of the variational EM algorithim. The \code{findLabels} function selects an 
#' optimal permutation of the labels to match a given ground truth.
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
#' @param training the ground truth theta
#' @param exhaustive a boolean for whether an exhaustive search should be performed. If false, a greedy algorithim is used instead
#' @return perm optimal permutation with respect to squared error loss
#' @return loss the sum of squared error loss of the optimal permutation weighted by relative frequency
#' @seealso permuteLabels
#' @examples
#' 
#' \dontshow{
#' set.seed(123)
#' Total <- 50 # 50 Individuals
#' J <- 3 # 2 different variables
#' dist <- c("multinomial", "bernoulli", "rank") # distributions of each variable
#' Rj <- c(100, 5, 1) # 100 repeated measures of the multinomial, 5 repeated measures of the bernoulli, 1 repeated measure of the rank
#' Nijr <- array(1, dim = c(Total, J, max(Rj))) #Nijr will always be 1 for multinomials and bernoulli's
#' K <- 4 # 4 sub-populations
#' alpha <- rep(.5, K) #hyperparameter for dirichlet distribution
#' Vj <- c(10, 1, 4) # Number of choices for each variable. Note the Bernoulli must have 1 choice
#' 
#' theta <- array(0, dim = c(J, K, max(Vj)))
#' # Parameters governing multinomial
#' theta[1,,] <- gtools::rdirichlet(K, rep(.3, Vj[1]))
#' #parameters governing bernoulli
#' theta[2,,] <- cbind(rbeta(K, 1,1), matrix(0, nrow = K, ncol = Vj[1]-1))
#' theta[3,,] <- cbind(gtools::rdirichlet(K, rep(.3, Vj[3])), matrix(0, nrow = K, ncol = Vj[1]-Vj[3]))
#' 
#' # generate group memberships
#' lambda <- gtools::rdirichlet(Total, rep(.2,K))
#' 
#' # Candidates selected for each observation. For Multinomial and Bernoulli, this is always 1
#' Nijr = array(0, dim = c(Total, J, max(Rj)))
#' Nijr[,1,c(1:Rj[1])] = 1
#' Nijr[,2,c(1:Rj[2])] = 1
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
#'       sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) # sub-population which governs the multinomial
#'       if(dist[j] =="multinomial") {
#'         obs[i,j,r,1] <- sample.int(Vj[j], size = 1, prob = theta[j,sub.pop,c(1:Vj[j])])-1 # must be in 0:(Vj[j]-1)
#'       }
#'       if(dist[j] == "rank") {
#'         obs[i,j,r,c(1:Nijr[i,j,r])] <- sample.int(Vj[j], size = Nijr[i,j,r], replace = FALSE, prob = theta[j,sub.pop,c(1:Vj[j])])-1  # must be in 0:(Vj[j]-1)
#'       }
#'       if(dist[j] == "bernoulli"){
#'         obs[i,j,r,1] <- rbinom(n = 1, size = 1, prob = theta[j,sub.pop,1])
#'       }
#'     }
#'   }
#' } 
#' }
#' 
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
        loss.i = sum(sweep((fitted.set[,perms[i,],]-training[,c(1:K),])^2, MARGIN = 2,weight, '*'))
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
#'Mixed Membership models are invariant to permutations of the sub-population labels and the ordering of the labels
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
  for(j in 1:out$J)
  {
    out$theta[j,,] = out$theta[j,perm,]
  }
  
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
  return(out)
}