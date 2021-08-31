
# These are a collection of convenience functions for use in numerical simulations
#
#    onesample.t.stat :  calculates the typical one-sample t statistic
#    rft.critical.t   :  calculates the RFT-corrected critical t threshold
#    rnorm1d          :  generates random Gaussian samples


onesample.t.stat <- function(y){
    m     <- colMeans( y )
    s     <- apply(y, 2, sd)
    n     <- nrow(y)
    t     <- m / ( s / sqrt(n) )
    return( t )
}


rft.critical.t <- function( r, alpha=0.05 ){
    # Calculate the critical t threshold using a Random Field Theory correction
    #    INPUTS
    #        r     :  (J,Q) array of residuals;  J = sample size, Q = number of domain nodes
    #        alpha :  Type I error rate
    #
    #    OUTPUTS
    #        tcrit :  critical t value
    #                 t values calculated from smooth 1D Gaussian data are
    #                 expected to reach tcrit with a probability of alpha
    #
    r   <- y - colMeans(y)                  # residuals
    df  <- nrow(y) - 1                      # degrees of freedom
    rn  <- r / sqrt( colSums( r^2 ) )       # normalized residuals
    rdn <- apply(rn, 1, diff)               # derivatives of normalized residuals
    L0  <- 1                                # 1 for unbroken fields
    L1  <- sum( sqrt( colSums( rdn^2 ) ) )  # Lipschitz-Killing curvature of the residuals
    eec <- function(u){                     # expected Euler characteristic
        a0 <- L0 * ( 1 - pt( u, df = df ) )
        a1 <- L1 * ( 1 + u^2 / df )^( ( 1 - df ) / 2 ) / ( 2 * pi )
        return( a0 + a1 )
    }
    tcrit <- uniroot( function(u) eec(u) - alpha, lower=1, upper=10 )$root  # threshold "u" that yields eec=alpha
    return( tcrit )
}


rnorm1d <- function( J=8, Q=101, nbasis=10, norder=4 ){
    library(fda)
    x       <- seq(0,1,length.out=Q)
    basis   <- create.bspline.basis(rangeval = c(0,1), nbasis=nbasis, norder=norder)
    coef    <- matrix( rnorm(J*nbasis), nbasis, J)
    fdobj   <- fd(coef, basis)
    y       <- eval.fd(x, fdobj)
    return( t(y) )
}

