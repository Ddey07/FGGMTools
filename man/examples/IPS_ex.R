require(igraph)
require(mvtnorm)
# Generate data y
p <- 10
Q <- qr.Q(qr(matrix(rnorm(p^2), p)))
Sigma <- crossprod(Q, Q*(5:1))
y = rmvnorm(1000, sigma=Sigma)
S = cov(y)
# Fix an adjacency matrix. Here we pick path graph of length p
A = as.matrix(as_adjacency_matrix(make_tree(p,1,"undirected")))
diag(A) = 1
## Solve
K = IPS(S,A)
# Check for the non-edges if they equal to zero in precision matrix
all(K[A==0]==0)
# Check if the edge-specific entries are equal to that of S
max(abs(solve(K)[A==1]  - S[A==1]))