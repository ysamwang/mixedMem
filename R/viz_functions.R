#' Mixed Membership Visualization
#' 
#' 
#' Plots \eqn{\theta}, the parameters which govern the distributions of variables in a mixed membership model. \eqn{\theta_{j,k}} governs the
#' distribution of variable j for sub-population k. Estimated distributions are shown with black bars. The ground truth theta (if available)
#' is shown in green dots. 
#'  
#' @param model the \code{mixedMemModel} object that will be plotted
#' @param truth optional argument if ground truth is available
#' @param main optional title argument
#' @param highlight matrix with 2 columns specifying any specific variable/sub-population combinations
#' to be highlighted in the plot  
#' @param varNames vector of length J specifying labels for each variable
#' @param groupNames vector of length K specifying labels for each sub-population
#' @export
vizTheta = function(model, truth = NULL, main = "Estimated Theta",
                    highlight = NULL, varNames = NULL, groupNames = NULL)
{
  par(oma = c(3,3,3,1))
  if(is.null(varNames))
  {
    varNames = paste("Var", c(1:model$J))
  }
  if(is.null(groupNames))
  {
    groupNames = paste("Group", c(1:model$K))
  }
  par(mfrow = c(min(model$J,10), model$K), mar = rep(.1,4))
  for(j in 1:model$J)
  {
    if(model$dist[j]=="multinomial"|model$dist[j]=="rank")
    {
      for(k in 1:model$K)
      {
        plot(model$theta[j,k,], type = "h", lwd = 2, col = "black", ylim = c(0,1), xlim = c(0, model$Vj[j]+.5),
             yaxt = "n", xaxt = "n")
        if(!is.null(truth))
        {
          points(c(1:model$Vj[j]), truth[j,k,], col = "blue", pch = 16)
        }
        if(k==1)
        {
          mtext(varNames[j], line = .2, side = 2, cex = .8)
        }
        if(j ==model$J| (j %%10)==0)
        {
          mtext(paste(groupNames[k], sep = " "), line = .2, side = 1, cex = 1-min(model$J,10)*.4)
        }
        if(!is.null(highlight))
        {
          for(r in 1:length(highlight[,1]))
          {
            if(j == highlight[r,1] & k == highlight[r,2])
            {
              box(which = "plot", lwd = 2, col = "red") 
            }
          }
        }
      }
    } else if (model$dist[j] == "bernoulli")
    {
      for(k in 1:model$K)
      {
        plot(model$theta[j,k,], type = "h", lwd = 2, col = "black", ylim = c(0,1), xlim = c(.5, 1.5),
             yaxt = "n", xaxt = "n")
        if(!is.null(truth))
        {
          points(c(1:model$Vj[j]), truth[j,k,], col = "green", pch = 19)
        }
        if(k==1)
        {
          mtext(varNames[j], line = .2, side = 2, cex = .8)
        }
        if(j ==model$J| (j %%10)==0)
        {
          mtext(paste(groupNames[k], sep = " "), line = .2, side = 1, cex = .8)
        }
        if(!is.null(highlight))
        {
          for(r in 1:length(highlight[,1]))
          {
            if(j == highlight[r,1] & k == highlight[r,2])
            {
              box(which = "plot", lwd = 2, col = "red") 
            }
          }
        }
      }
    } 
    if((j %%10) == 0)
    {
      title(main = main, outer = T, cex = 1.2)
    }
    
  }
  title(main = main, outer = T, cex = 1.2)
}
