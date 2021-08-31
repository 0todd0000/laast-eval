
# This script calculates the false positive rate for
#    a one-sample LAAST implementation: "mylaast.onesample"
 

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics




onesample.t.stat <- function(y){
    m     <- colMeans( y )
    s     <- apply(y, 2, sd)
    n     <- nrow(y)
    t     <- m / ( s / sqrt(n) )
    return( t )
}


rft.critical.t <- function( y, alpha=0.05 ){
    # Calculate the critical t threshold using a Random Field Theory correction
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




# (0) Check single sample generation
J      <- 8
Q      <- 101
nbasis <- 10
norder <- 3
y      <- rnorm1d(J, Q, nbasis, norder)
graphics.off()
matplot( t(y), type='l', lty=1, lwd=0.5, col='blue')



# (1) Source R files:
dirREPO  <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirR     <- file.path( dirREPO, 'R' )     # path to the R directory in this repository
# source the required laast-eval functions:
fnameR1  <- file.path(dirR, 'laast-eval', 'laast.R')
fnameR2  <- file.path(dirR, 'laast-eval', 'loess.R')
source( fnameR1 )
source( fnameR2 )





# (2) Simulate:

# set main simulation parameters:
set.seed(2)
J         <- 8   # sample size
Q         <- 100 # number of domain nodes
nbasis    <- 10  # number of basis functions (higher values yield rougher fields)
norder    <- 4   # bspline basis function order

# initialize simulation counts:
n         <- 0   # number of datasets generated
nfp.rft   <- 0   # number of false positives (RFT)
nfp.laast <- 0   # number of false positives (LAAST)

# run simulation:
set.seed(2)
for (i in 1:1000){
    n    <- n + 1 # total number of datasets
    y    <- rnorm1d(J, Q, nbasis, norder)
    
    # RFT inference (on the typical, least-squares, one-sample t statistic):
    t    <- onesample.t.stat(y)
    tc   <- rft.critical.t(y, alpha=0.05)  # critical t value
    if ( max(t) > tc ){
        nfp.rft <- nfp.rft + 1
    }
    
    # LAAST:
    results  <- mylaast.onesample(y, binSize=2, span=NULL, loess=T)
    pc       <- results[[1]]   # critical p value
    p        <- results[[2]]$pvals
    if ( min(p) < pc ){
        nfp.laast <- nfp.laast + 1
    }

    # report results
    msg <- sprintf("Iteration: %d, FPR (RFT): %.3f, FPR (LAAST): %.3f", n, nfp.rft/n, nfp.laast/n)
    print( msg )
}



