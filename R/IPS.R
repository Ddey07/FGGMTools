#' @title Get a graph-constrained MLE of the covariance matrix under Gaussian Graphical Model 
#'
#' @description Maximum likelihood estimate of the covariance matrix under the constraint of a Gaussian Graphical Model using Iterative Proportional Scaling algorithm. 
#'
#' @param S The sample covariance matrix (p by p), where p is the number of variables. 
#' @param A An adjacency matrix corresponding to the underlyning graphical model (p by p). Must come from an undirected graph. 
#' @param eps Error tolerance threshold (in likelihood) for the algorithm to converge
#' @return \code{IPS} returns the estimated precision matrix under the graphical model. 
#' @author Debangan Dey, Sudipto Banerjee, Martin Lindquist and Abhirup Datta
#' @references Dempster AP, Covariance selection - Biometrics, 1972 - JSTOR
#' @references DR Musgrove, J Hughes, LE Eberly, Hierarchical copula regression models for areal data- Spatial Statistics, 2016 - Elsevier
#' @references Ping-Feng Xu, Jianhua Guo and Man-Lai Tang, A Localized Implementation of the Iterative Proportional Scaling Procedure for Gaussian Graphical Models. Journal of Computational and Graphical Statistics.Vol. 24, No. 1 (MARCH 2015), pp. 205-229. 
#' @references Dey D., Banerjee S., Lindquist M., and Datta A., Graph-constrained Analysis for Multivariate Functional Data. Available at arXiv.org
#' @details
#' This code uses slightly modified functions from the testing version of copCS package available at: \url{https://github.com/donaldmusgrove/copCS/blob/master/R/ips.R} by Dr. D.R. Musgrove.
#' @export
#' @import stats
#' @example man/examples/IPS_ex.R

IPS <- function(S, A, eps=1e-8){
  n <- nrow(A)
  diag(A) <- 0
  T <- S
  
  W <- clique_part(A)
  Q <- IPS2(S=S, n=n, T=T, W=W, maxit=100, eps=eps)$Q
  
  return(Q)   
}
