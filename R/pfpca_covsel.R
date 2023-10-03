
#' @title Graph-constrained covariance estimation of multivariate function data 
#'
#' @description Estimates the graph-constrained Karhunen-Loeve expansion for a partially separable multivariate Gaussian process.
#' @param y  list of length p containing densely observed multivariate (p-dimensional) functional data . \code{y[[j]]} is an nxm matrix of functional data for n subjects observed on a grid of length m
#' @param t  (optional) grid on which functional data is observed, defaults to seq(0, 1, m) where \code{m = dim(data[[1]])[2]}
#' @param A  An adjacency matrix corresponding to the underlyning graphical model (p by p). Must come from an undirected graph. 
#' @param FVE Desired proportion of variance explained in the functional principal components analysis
#' @return A list of length L (number of princtipal components explaining desired FVE), each with two variables:
#' \describe{
#'   \item{\code{phi}}{A row denotes the value of l-th basis function evaluated at a grid of length m}
#'   \item{\code{Sigma_l}}{Estimated covariance matrix for the l-th basis coefficients under graphical constraint}
#' }
#' @note From the above list, the endowed covariance function can be obtained as \eqn{\sum_{l=1}^{L} \Sigma_l \phi_l \phi_l^T}. 
#' @author Debangan Dey, Sudipto Banerjee, Martin Lindquist and Abhirup Datta
#' @references Dey D., Banerjee S., Lindquist M., and Datta A., Graph-constrained Analysis for Multivariate Functional Data. Available at arXiv.org
#' @references Zapata, J., Oh, S., and Petersen, A. (2019) Functional Graphical Models for Partially Separable Gaussian Processes. Available at arXiv.org
#' @details
#' This code uses slightly modified functions \code{pfpca} and some helper functions from the fgm package available at: \url{https://github.com/javzapata/fgm} by Zapata, J., Oh, S., and Petersen, A.
#' @export
#' @import stats
#' @example man/examples/pfpca_covsel_ex.R

pfpca_covsel=function(y, t=seq(0, 1, length.out=dim(y[[1]])[2]), FVE ,A){
  
  # run pfpca 
  pf = pfpca(y,t, FVE)
  L = pf$L
  phimat = pf$phi
  psepcov = lapply(1:L,function(x){Sigma_l = cov(t(pf$theta[[x]]))
  # stabilize condition number of sigma_l 
  Sigma_l = as.matrix(nearPD(Sigma_l, posd.tol = 0.001)$mat)
  Sigma_l1 = IPS(S=Sigma_l, A= A, eps = 1e-4)
  Sigma_l1 = solve(Sigma_l1)
  list(Sigma_l = Sigma_l1, phi = phimat[x,])
  #kronecker(Sigma_l1, phimat[x,] %*% t(phimat[x,]))
  })
  return(psepcov)
}


pfpca = function(y, t=seq(0, 1, length.out=dim(y[[1]])[2]), FVE){

  # checking inputs
  if (missing(y)) stop('y is missing')
  if (!is.list(y)) stop('y must be a list of matrices')
  if (length(t) != dim(y[[1]])[2]) stop('length of time grid does not match dim(y[[1]])[2])')

  # Compute autocovariance surface of each component separately (smoothed first)
  p = length(y)
  optns = list(dataType = 'Dense', error = FALSE)
  mu = t(sapply(y, function(yy) .GetMeanDense(yy, t, optns)$mu))
  m = length(t); n = dim(y[[1]])[1]
  S = array(0, dim = c(m, m, p))
  for(j in 1:p){
    S[,,j] = .GetCovDense(y[[j]], mu[j,], optns)$smoothCov
  }
  # Add together and compute eigendecomposition

  optns$FVEthreshold = FVE
  optns$maxK = n*p
  C = rowSums(S, dims = 2)
  E = .GetEigenAnalysisResults(C, t, optns)

  # Compute vector scores

  y2 = lapply(1:n, function(i){ # Rearrange data by subject
    t(sapply(1:p, function(j) y[[j]][i,]))
  })

  L = dim(E$phi)[2]

  #if(L > dim(E$phi)[2]){
  #  #warning('L is too large, changing to maximum value')
  #  L = dim(E$phi)[2]
  #}

  theta = lapply(1:L, function(l){
    sapply(1:n, function(i){

      sapply(1:p, function(j) fdapace::trapzRcpp(t, (y2[[i]][j,] - mu[j,])*(E$phi[,l])))
    })
  })

  # Compute FVE
  FVE = E$cumFVE[L]

  return(list(phi = t(E$phi[,1:L]), theta = theta, FVE= E$cumFVE[1:L], L=L))
  #FVE = FVE,
  #cumFVE=E$cumFVE, L=L))

}

.GetEigenAnalysisResults <- function(smoothCov, regGrid, optns, muWork = NULL) {
  # this function is based on <https://github.com/functionaldata/tPACE/blob/master/R/GetEigenAnalysisResults.R>
  maxK <- optns$maxK
  FVEthreshold <- optns$FVEthreshold
  verbose <- optns$verbose

  gridSize <- regGrid[2] - regGrid[1]
  numGrids <- nrow(smoothCov)

  eig <- eigen(smoothCov)

  positiveInd <- eig[['values']] >= 0
  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }
  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]

  if (maxK < length(d)) {
    if (optns[['verbose']]) {
      message(sprintf("At most %d number of PC can be selected, thresholded by `maxK` = %d. \n", length(d), maxK))
    }

    d <- d[1:maxK]
    eigenV <- eigenV[, 1:maxK, drop=FALSE]
  }

  # thresholding for corresponding FVE option
  #(not before to avoid not being able to reach the FVEthreshold when pos eigenvalues > maxk)
  # i.e. default FVE 0.9999 outputs all components remained here.
  FVE <- cumsum(d) / sum(d) * 100  # cumulative FVE for all available eigenvalues from fitted cov
  no_opt <- min(which(FVE >= FVEthreshold * 100)) # final number of component chosen based on FVE

  # normalization
  if (is.null(muWork)) {
    muWork = 1:dim(eigenV)[1]
  }

  phi <- apply(eigenV, 2, function(x) {
    x <- x / sqrt(fdapace::trapzRcpp(regGrid, x^2))
    #x <- x / sqrt(pracma:::trapz(regGrid, x^2))
    if ( 0 <= sum(x*muWork) )
      return(x)
    else
      return(-x)
  })
  lambda <- gridSize * d;

  fittedCov <- phi %*% diag(x=lambda, nrow = length(lambda)) %*% t(phi)

  return(list(lambda = lambda[1:no_opt], phi = phi[,1:no_opt, drop=FALSE],
              cumFVE = FVE, kChoosen=no_opt, fittedCov=fittedCov))
}


.GetMeanDense <- function(ymat, obsGrid, optns){
  # This function is based on <https://github.com/functionaldata/tPACE/blob/master/R/GetMeanDense.R>
  # Check optns
  if(!(optns$dataType %in% c('Dense', 'DenseWithMV'))){
    stop('Cross sectional mean is only applicable for option: dataType = "Dense" or "DenseWithMV"!')
  }

  if ( is.null(optns$userMu) ){
    mu = colMeans(ymat, na.rm = TRUE) # use non-missing data only
  } else {
    mu = spline(optns$userMu$t, optns$userMu$mu, xout= obsGrid)$y;
  }

  if(any(is.na(mu))){
    stop('The cross sectional mean is appears to have NaN! Consider setting your dataType to \'Sparse\' manually')
  }

  ret = list('mu' = mu, 'muDense' = NULL, 'mu_bw' = NULL)
  class(ret) = "SMC"
  return(ret)
}




.GetCovDense <- function(ymat, mu, optns) {
  # This function is based on <https://github.com/functionaldata/tPACE/blob/master/R/GetCovDense.R>
  if(!(optns$dataType %in% c('Dense', 'DenseWithMV'))){
    stop('Sample Covariance is only applicable for option: dataType = "Dense" or "DenseWithMV"!')
  }
  # if( optns$muCovEstMethod == 'cross-sectional' ){
  n = nrow(ymat)
  m = ncol(ymat)
  if( !is.null(optns$userMu) ){
    ymat = ymat - matrix(rep(times= nrow(ymat), mu), ncol= ncol(ymat), byrow=TRUE)
    K = matrix( rep(0,m^2), m)
    for( i in (1:m)){
      for( j in (1:m)){
        XcNaNindx = which(is.na(ymat[,i]));
        YcNaNindx = which(is.na(ymat[,j]));
        NaNrows = union(XcNaNindx, YcNaNindx);
        # Find inner product of the columns with no NaN values
        indx = setdiff( 1:n, NaNrows)
        K[i,j] =  sum(ymat[indx,i] * ymat[indx,j]) * (1/(n-1-length(NaNrows)));
      }
    }
  } else {
    K = cov(ymat, use = 'pairwise.complete.obs') # sample variance using non-missing data
  }
  K = 0.5 * (K + t(K)) # ensure that K is symmetric

  if (any(is.na(K))) {
    stop("Data is too sparse to be considered DenseWithMV. Remove sparse observations or specify dataType='Sparse' for FPCA")
  }

  if (optns[['error']] == TRUE) {
    # 2nd order difference method for finding sigma2
    if (!is.null(optns[['userSigma2']])) {
      sigma2 <- optns[['userSigma2']]
    } else {
      #browser()
      ord <- 2
      sigma2 <- mean(diff(t(ymat), differences=ord)^2, na.rm=TRUE) /
        choose(2 * ord, ord)
      diag(K) <- diag(K) - sigma2
    }
  } else {
    sigma2 <- NULL
  }

  ret = list('rawCov' = NULL, 'smoothCov' = K, 'bwCov' = NULL,
             'sigma2' = sigma2, outGrid = NULL)
  class(ret) = "SmoothCov"
  return(ret)
  # } else {
  # stop('optns$muCovEstMethod is unknown!\n')
  # }
}
