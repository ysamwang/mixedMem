# Internal function to check consistency for mixedMemModelVI class
checkModelVarInf <- function(model)
{
  Total <- model[[1]]
  J <- model[[2]]
  Rj <- model[[3]]
  Nijr <- model[[4]]
  K <- model[[5]]
  Vj <- model[[6]]
  alpha <- model[[7]]
  theta <- model[[8]]
  phi <- model[[9]]
  delta <- model[[10]]
  dist <- tolower(model[[11]])
  obs <- model[[12]]

  if (class(model) != "mixedMemModelVarInf") {
    stop("Input must be: mixedMemModelVarInf")
  }

  #Check Total
  if (Total < 1 || Total != round(Total)) {
    stop("Input must be positive Integer: ", "Total")
  }
  if (length(Total) > 1) {
    stop("Input of incorrect dimensions: ", "Total", " must be of dimension ", "1")
  }
  
  # Check J
  if (J < 1 || J != round(J)) {
    stop("Input must be positive Integer: ", "J")
  }
  if (length(J) > 1) {
    stop("Input of incorrect dimensions: ", "J", " must be of dimension ", "1")
  }
  
  
  # Check Rj
  if (Rj < 1 || Rj != round(Rj)) {
    stop("Input must be positive Integer: ", "Rj")
  }
  
  if (length(Rj) != J | !is.null((dim(Rj)))) {
    stop("Input of incorrect dimensions: ", "Rj", " must be of dimension ", "J")
  }
  maxR = max(Rj)
  
  # Check Nijr
  if ( any(Nijr < 1) | any(Nijr != round(Nijr))) {
    stop("Input must be positive Integer: ", "Nijr")
  }
  if (any(dim(Nijr) != c(Total,J,maxR))) {
    stop("Input of incorrect dimensions: ", "Nijr", " must be of dimension ", "{Total,J,max(Rj)}")
  }
  maxN = max(Nijr)
  
  # Check K
  if (K < 1 | K != round(K)) {
    stop("Input must be positive Integer: ", "K")
  }
  if (length(K) > 1) {
    stop("Input of incorrect dimensions: ", "K", " must be of dimension ", "1")
  }
  
  # check dist
  if (!all(dist %in% c("bernoulli", "multinomial", "rank"))) {
    stop("dist is invalid! Valid choices are bernoulli, multinomial, or rank")
  }
  if (length(dist) != J) {
    stop("dist must be of length J")
  }
  
  # Check Vj
  if (Vj < 1 || Vj != round(Vj)) {
    stop("Input must be positive Integer: ", "Vj")
  }
  if (length(Vj) != J) {
    stop("Input of incorrect dimensions: ", "Vj", " must be of dimension ", "J")}
  if (Vj > 1 && dist == "bernoulli") {
    stop("V must be 1 for bernoulli variables")
  }
  
  for (j in 1:J) {
    if (any(Nijr[,j,] > Vj[j])) {
      stop(paste("For variable ",j, " N_ijr must be less than V_j"))
    }
  }

  maxV = max(Vj)
  
  # Check alpha
  if (any(alpha < 0) ) {
    stop("Input must be positive: ", "alpha")
  }
  if (length(alpha) != K | !is.null(dim(alpha)) ) {
    stop("Input of incorrect dimensions: ", "alpha", " must be of dimension ", "K")
  }
  
  #check theta
  if (any(theta < 0) ) {
    stop("Input must be positive: ", "theta")
  }
  if (any(dim(theta) != c(J,K,maxV))) {
    stop("Input of incorrect dimensions: ", "theta", " must be of dimension ", "{1,J,K,max(V)}")}
  
  for (j in 1:J) {
    for (k in 1:K) {
      if (abs(sum(theta[j,k,]) - 1) > 1e-8 & dist[j] != "bernoulli") {
        stop("Distribution must sum to 1 for theta for Variable ", j)
      }
      if (any(theta[j,k,] > 1 )) {
        stop("Distribution must sum to 1 for theta for Variable ", j)
      }
      if (any(theta[j,k, 1:Vj[j]] < 1e-13) ) {
        stop("Theta must be bounded away from 0 (for numerical stability)")
      }
    }
  }
  
  
  #check phi
  if (any(phi < 0) ) {
    stop("Input must be positive: ", "phi")
  }
  if (any(dim(phi) != c(Total,K))) {
    stop("Input of incorrect dimensions: ", "phi", " must be of dimension ", "{Total,K}")
  }
  
  
  #check delta
  if (any(delta < 0) ) {
    stop("Input must be positive: ", "delta")
  }
  if (any(dim(delta) != c(Total,J,maxR,maxN,K))) {
    stop("Input of incorrect dimensions: ", "delta", " must be of dimension ", "{Total,J,max(Rj),max(Nijr),K}")
  }
  
  for (i in 1:Total) {
    for (j in 1:J) {
      for (r in 1:Rj[j]) {
        for (n in 1:Nijr[i,j,r]) {
          if (abs(sum(delta[i,j,r,n,]) - 1) > 1e-10) {
            stop("Delta must sum to 1 for delta[",paste(i,j,r,n,sep = ","),", ]")
          }
        }
      }
    }
  }
  
  #check obs
  if (obs != round(obs) || obs < 0 )   {
    stop("obs must be non-negative integers")
  }
  if (any(dim(obs) != c(Total, J, maxR, maxN))) {
    stop("Input of incorrect dimensions: ", "Obs", " must be of dimension ", "{Total,J,max(Rj),max(Nijr)}")
  }
}



# Internal function to check mixedMemModelMCMC object
checkModelMCMC <- function(model){ 
    Total <- model[[1]]
    J <- model[[2]]
    dist <- tolower(model[[3]])
    
    Rj <- model[[4]]
    Vj <- model[[5]]
    obs <- model[[6]]
    
    K <- model[[7]]
    
    theta <- model[[8]]
    alpha <- model[[9]]
    ksi <- model[[10]]
    
    lambda <- model[[11]]
    Z <- model[[12]]
    
    tau <- model[[13]]
    beta <- model[[14]]
    gamma <- model[[15]]
    extended <- model[[16]]
    if (extended) {
      P <- model[[17]]
      S <-  model[[18]]
      fixedObs <- model[[19]]
    }
    
    
    #Check Total
    if (Total < 1 || Total != round(Total)) {
      stop("Input must be positive Integer: ", "Total")
    }
    if (length(Total) > 1) {
      stop("Input of incorrect dimensions: ", "Total", " must be of dimension ", "1")
    }
    
    # Check J
    if (J < 1 || J != round(J)) {
      stop("Input must be positive Integer: ", "J")
    }
    if (length(J)>1) {
      stop("Input of incorrect dimensions: ", "J", " must be of dimension ", "1")
    }
    
    
    # Check Rj
    if (Rj < 1 || Rj != round(Rj)) {
      stop("Input must be positive Integer: ", "Rj")
    }
    if (length(Rj) != J|!is.null((dim(Rj)))) {
      stop("Input of incorrect dimensions: ", "Rj", " must be of dimension ", "J")
    }
    maxR = max(Rj)
    
    
    # Check K
    if (K < 1 | K != round(K)) {
      stop("Input must be positive Integer: ", "K")
    }
    if (length(K) > 1) {
      stop("Input of incorrect dimensions: ", "K", " must be of dimension ", "1")
    }
    
    # check dist
    if (!all(dist %in% c("bernoulli", "multinomial"))) {
      stop("dist is invalid! Valid choices are bernoulli or multinomial")
    }
    if (length(dist) != J) {
      stop("dist must be of length J")
    }
    
    # Check Vj
    if (Vj < 1 || Vj != round(Vj)) {
      stop("Input must be positive Integer: ", "Vj")
    }
    if (length(Vj) !=J) {
      stop("Input of incorrect dimensions: ", "Vj", " must be of dimension ", "J")}
    if (Vj > 1 && dist=="bernoulli") {
      stop("V must be 1 for bernoulli variables")
    }
    maxV = max(Vj)
    
    # Check alpha
    if (alpha < 0) {
      stop("Input must be positive: ", "alpha")
    }
    
    if( length(alpha) != 1){
      stop("Alpha must be scalar valued")
    }
    
    #check theta
    if (any(theta < 0) ) {
      stop("Input must be positive: ", "theta")
    }
    
    if (any(dim(theta)!=c(J,K,maxV))) {
      stop("Input of incorrect dimensions: ", "theta", " must be of dimension ", "{1,J,K,max(V)}")}
    
    for(j in 1:J)
    {
      for(k in 1:K)
      {
        if (abs(sum(theta[j,k,])-1) > 1e-8 & dist[j] !="bernoulli") {
          stop("Distribution must sum to 1 for theta for Variable ", j)
        }
        if (any (theta[j,k,] > 1 )) {
          stop("Distribution must sum to 1 for theta for Variable ", j)
        }
      }
    }
    
    
    # check phi
    if (beta < 0) {
      stop("Input must be positive: ", "beta")
    }
    
    # check delta
    if (gamma < 0) {
      stop("Input must be positive: ", "gamma")
    }
    
    
    
    #check obs
    if (obs != round(obs) || obs < 0)   {
      stop("obs must be non-negative integers")
    }
    if (any(dim(obs)!= c(Total, J, maxR))) {
      stop("Input of incorrect dimensions: ", "Obs", " must be of dimension ", "{Total,J,max(Rj)}")
    }
    
    # If extended model
    if (extended) {
      # check fixed obs
      if(S < 1){
        stop("Extended model specified: S must be greater than 1")
      }
      
      if (!is.null(fixedObs) ){
        if (any(dim(fixedObs) != c(S, J, maxR))) {
          stop("Extended model specified:fixedObs must be of dimension (S, J, max(Rj))")
        }
      } else {
        stop("Extended model specified: fixedObs must be specified")
      }
      
      if(is.null(P)) {
        stop("Extended model specified: P must be specified")
      } else if(length(P) != (S + 1)) {
        stop("Extended model specified: P must be of length S + 1")
      }
    }
    
}
