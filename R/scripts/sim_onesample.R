
# This script calculates the false positive rate for
#    a one-sample LAAST implementation: "mylaast.onesample"
 

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics




# The GaussDensitySumNoise function below is a
#   random (Gaussian) function generator that appears in
#   Fabian Telschow's "SCBfda" repository on GitHub:
#   https://github.com/ftelschow/SCBfda
#
#   The code and comments below were copy-and-pasted from:
#      https://github.com/ftelschow/SCBfda/blob/master/R/RandomFieldGeneration.R

#' Creates sample paths from a 1D field generated as a random sum of Gaussian densities with different means and variances and random coefficients.
#' @param N Integer amount of realisations of the random field to be generated.
#' @param x Vector locations at which the random field is evaluated.
#' @param sigma Function computing the pointwise variance of the field. Default value is unit variance everywhere.
#' @param randNumber Function generating a vector of random numbers with mean zero and variance 1. Default is rnorm().
#' @param anchorPoints Vector containing the locations of the mean of the Gaussian densities in the random sum. Default value are 15 equidistant knots of the interval specifying the domain.
#' @param anchorSd Vector containing the standard deviations of the stand
#' @return Matrix containing the realisations of the random field as columns. The rows are the locations.
#' @export
GaussDensitySumNoise <- function(N, x = seq(0,1,length.out=100),
                      sigma = function(x){ rep(1, length(x)) }, randNumber = rnorm, anchorPoints = NULL, anchorSd = NULL){
  if(is.null(anchorPoints)){
    anchorPoints <- seq(from=min(x), to=max(x), length.out=15)
  }
  nAnchorPoints = length(anchorPoints)
  if(is.null(anchorSd)){
    anchorSd <- rep(diff(range(x))/nAnchorPoints, nAnchorPoints)
  }
  f <- sapply(1:length(anchorPoints),
              function(pt) dnorm(x, mean=anchorPoints[pt],
                                 sd=anchorSd[pt]))
  fSqSum <- apply(f^2, 1, sum)
  fNorm <- f / sqrt(fSqSum)
  y <- (fNorm %*% matrix(randNumber(nAnchorPoints * N), nAnchorPoints, N)) * sigma(x)
  return( t(y) )
}


onesample.t.stat <- function(y){
    m     <- colMeans( y )
    s     <- apply(y, 2, sd)
    n     <- nrow(y)
    t     <- m / ( s / sqrt(n) )
    return( t )
}


rft.critical.t <- function( y, alpha=0.05 ){
    # critical t threshold (Random Field Theory correction)
    r   <- y - colMeans(y)  # residuals
    df  <- nrow(y) - 1      # degrees of freedom
    rdn <- diff( r / sqrt( colSums( r^2 ) ) ) # derivatives of normalized residuals
    L1  <- mean( colSums( rdn^2 ) )  # Lipschitz-Killing curvature of the residuals
    L0  <- 1                # first "resel count";  1 for unbroken fields
    eec <- function(u){     # expected Euler characteristic
        a0 <- L0 * ( 1 - pt( u, df = df ) )
        a1 <- L1 * ( 1 + u^2 / df )^( ( 1 - df ) / 2 ) / ( 2 * pi )
        return( a0 + a1 )
    }
    tcrit <- uniroot( function(u) eec(u) - alpha, lower=1, upper=10 )$root
    return( tcrit )
}


randn1d <- function( J=8, Q=101, fwhm=20 ){
    library(mmand)
    y      <- matrix( rnorm(J*Q), nrow=J, ncol=Q)
    sigma  <- fwhm / sqrt( 8 * log(2) )
    y      <- t( gaussianSmooth( t(y), sigma) )
    ### define Gaussian kernel
    t          = c(  (-0.5*(Q-1)) : (0.5*(Q-1))  )
    gf         = exp(-(t^2) / (2*sigma^2))
    gf         = gf / sum(gf)
    ### expected variance for this kernel
    AG         = fft(gf)
    Pag        = AG * Conj(AG)  #power of the noise
    COV        = Re( fft(Pag, inverse=T) / Q )
    svar       = COV[1]
    scale      = sqrt( 1 / svar )
    return ( y * scale )
}



# library(fda.usc)


# J <- 8
# Q <- 100
# FWHM <- 20
# y <- randn1d(J, Q, FWHM)




# graphics.off()
# matplot( t(y), type='l', lty=1, lwd=0.5, col='blue')
# matplot( t(randn1d(8,101,20)) , type='l', lty=1, lwd=0.5, col='blue')



# # set.seed(3)
# n   <- 8
# y   <- GaussDensitySumNoise(n, x = seq(0,1,length.out=100))
# r   <- y - colMeans(y)  # residuals
# rdn <- diff( r / sqrt( colSums( r^2 ) ) ) # derivatives of normalized residuals
# lkc <- mean( colSums( rdn^2 ) )  # Lipschitz-Killing curvature of the residuals





# assemble directories:
dirREPO  <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirR     <- file.path( dirREPO, 'R' )     # path to the R directory in this repository


# source the required laast-eval functions:
fnameR1  <- file.path(dirR, 'laast-eval', 'laast.R')
fnameR2  <- file.path(dirR, 'laast-eval', 'loess.R')
source( fnameR1 )
source( fnameR2 )







# # simulate:
# set.seed(2)
# J <- 8  # sample size
# Q <- 100
# FWHM <- 5
# n   <- 0  # number of datasets
# nfp <- 0  # number of false positives
#
# tc <- 2.0605586830770948
# for (i in 1:1000){
#     # y        <- GaussDensitySumNoise(ss)
#     y <- randn1d(J, Q, FWHM)
#
#     # s     <- apply(y, 2, sd)
#     t        <- onesample.t.stat(y)
#     tc       <- rft.critical.t(y, alpha=0.05)
#
#     n     <- n + 1
#     fp    <- ( max( t ) > tc)
#     # print( max())
#     if (fp){
#         nfp <- nfp + 1
#     }
#
#
#
#     # results  <- mylaast.onesample(y, binSize=2, span=NULL, loess=T)
#     # pcrit    <- results[[1]]
#     # df           <- results[[2]]
#     # n     <- n + 1
#     # fp           <- min( df$pvals ) < pcrit
#     # if (fp){
#     #     nfp <- nfp + 1
#     # }
#     # print( c('Number of datasets: ', n, nfp / n) )
#     msg <- sprintf("Number of datasets: %d, False positive rate: %.3f", n, nfp/n)
#     print( msg )
# }



