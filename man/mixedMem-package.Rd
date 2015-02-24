\name{mixedMem-package}
\alias{mixedMem-package}
\alias{mixedMem}
\docType{package}

\title{
Tools for fitting discrete multivariate mixed membership models
}

\description{
The \code{mixedMem} package contains tools for fitting and interpreting discrete multivariate mixed membership models following the general framework outlined in Erosheva 2004. In a mixed membership models, individuals can
belong to multiple groups instead of only a single group (Airoldi 2014). This extension allows for a richer description of heterogeneous populations and has been applied in a wide variety of contexts including: text data (Blei 2003), genotype sequences (Pritchard 2000), ranked votes (Gormley 2009) and survey data (Erosheva 2007).
}

\details{
Mixed membership model objects can be created using the \code{mixedMemModel}
constructor function. This function checks the internal consistency of the data/parameters
and returns an object suitable for use by the \code{mmVarFit} function. The 
\code{mmVarFit} function is the main function in the package. It utilizes a variational EM algorithim to fit an approximate posterior distribution for the latent variables and select pseudo-MLE estimates for the hyperparameters. A step-by-step guide to using the package is detailed in the package vignette. 

The package supports multivariate models (with or without repeated measurements) where each variable can be of a different type. Currently supported data types include: Bernoulli, rank (Plackett-Luce) and multinomial. We assume the following generative model for each mixed membership model:

For each individual i = 1,\ldots Total: Draw \eqn{\lambda_i} from a Dirichlet(\eqn{\alpha}). \eqn{\lambda_i} is a vector of length K which indicates the degree of membership for individual i in each of the K sub-populations.
  \itemize{
  \item{For each variable j = 1, \ldots, J, each of replicate r = 1, \dots, \eqn{R_j} and each ranking level n = 1\ldots, \eqn{N_{i,j,r}}-
  
  Draw \eqn{Z_{i,j,r,n}} from a multinomial(1, \eqn{\lambda_i}). \eqn{Z_{i,j,r,n}} determines the sub-population which governs the response for observation \eqn{X_{i,j,r,n}}. This is sometimes referred to as the context vector because it determines the context from whicht the individual responds.}
  \item{For each variable j = 1, \ldots, J, each of replicate r = 1, \dots, R_j and each ranking level n = 1\ldots, \eqn{N_{i,j,r}}-
  
  Draw \eqn{X_{i,j,r,n}} from the distribution parameterized by \eqn{\theta_{j,Z_{i,j,r,n}}}. \eqn{\theta} is the set of parameters which govern the observations for each sub-population. If variable j is a multinomial or rank distribution with \eqn{V_j} categories/candidates, \eqn{\theta_{j,k}} is a vector of length \eqn{V_j} which parameterizes the responses to variable j for sub-population k. If variable j is a Bernoulli random variable, then \eqn{\theta_{j,k}} is a value which determines the probability of success. 
}
}
}
\author{
Sam Wang <ysamwang@uw.edu>, Elena Erosheva <Elena.A.Erosheva@uw.edu>
}

\references{
Airoldi, E. M., Blei, D., Erosheva, E. A., & Fienberg, S. E. (2014). Handbook of Mixed Membership Models and Their Applications. CRC Press. Chicago  

Blei, David M.; Ng, Andrew Y.; Jordan, Michael I. (2003). Latent Dirichlet Allocation. Journal of Machine Learning Research, 3, 993-1022. \url{http://www.cs.columbia.edu/~blei/papers/BleiNgJordan2003.pdf}

Erosheva, Elena A.; Fienberg, Stephen E.; Joutard, Cyrille (2007). Describing disability through individual-level mixture models for multivariate binary data. The Annals of Applied Statistics 1 (2007), no. 2, 502--537. doi:10.1214/07-AOAS126. \url{http://projecteuclid.org/euclid.aoas/1196438029}.

Erosheva, Elena A.; Fienberg, Stephen E.; Lafferty, John. (2004). Mixed-membership models of scientific publications".
PNAS, 101 (suppl 1), 5220-5227. doi:10.1073/pnas.0307760101. \url{http://www.pnas.org/content/101/suppl_1/5220.full}.

Gormley, Isobel C.; Murphy, Thomas B. (2009). A grade of membership model for rank data. Bayesian Analysis, 4, 265 - 296. DOI:10.1214/09-BA410. \url{http://ba.stat.cmu.edu/journal/2009/vol04/issue02/gormley.pdf}

National Election Studies, 1983 Pilot Election Study. Ann Arbor, MI: University of Michigan, Center for Political Studies, 1999

Pritchard, Jonathan K., Matthew Stephens, and Peter Donnelly. Inference of population structure using multilocus genotype data. Genetics 155.2 (2000): 945-959.

}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{mixed membership, grade of membership}
