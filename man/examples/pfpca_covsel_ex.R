## Variables
# Omega - list of precision matrices, one per eigenfunction
# Sigma - list of covariance matrices, one per eigenfunction
# theta - list of functional  principal component scores
# phi - list of eigenfunctions densely observed on a time grid
# y - list containing densely observed multivariate (p-dimensional) functional data

library(mvtnorm)
library(fda)
library(igraph)


## Generate data y
 source(system.file("exec", "getOmegaSigma.R", package = "fgm"))
 theta = lapply(1:nbasis, function(b) t(rmvnorm(n = 100, sigma = Sigma[[b]])))
 theta.reshaped = lapply( 1:p, function(j){
     t(sapply(1:nbasis, function(i) theta[[i]][j,]))
 })
 phi.basis=create.fourier.basis(rangeval=c(0,1), nbasis=21, period=1)
 t = seq(0, 1, length.out = time.grid.length)
 chosen.basis = c(2, 3, 6, 7, 10, 11, 16, 17, 20, 21)
 phi = t(predict(phi.basis, t))[chosen.basis,]
 y = lapply(theta.reshaped, function(th) t(th)%*%phi)

# Fix an adjacency matrix for the graph between the variables
A = as_adjacency_matrix(make_tree(length(y),1,"undirected"))
# Get graph-constrained estimate of covariance function for the process
pf_ips = pfpca_covsel(y,A=A,FVE=0.8)
