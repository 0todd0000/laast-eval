
# This script calculates the false positive rate for
#    a one-sample LAAST implementation: "mylaast.onesample"
 

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics





# (0) Source R files:
dirREPO  <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirR     <- file.path( dirREPO, 'R' )     # path to the R directory in this repository
# source the required functions:
fnameR1  <- file.path(dirR, 'laast-eval', 'sim.R')
source( fnameR1 )




# (1) Check single sample generation
J      <- 8
Q      <- 101
nbasis <- 10
norder <- 3
y      <- rnorm1d(J, Q, nbasis, norder)
graphics.off()
matplot( t(y), type='l', lty=1, lwd=0.5, col='blue')





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

# run simulation:
set.seed(2)
for (i in 1:1000){
    n    <- n + 1 # total number of datasets
    y    <- rnorm1d(J, Q, nbasis, norder)
    
    # RFT inference (on the typical, least-squares, one-sample t statistic):
    t    <- onesample.t.stat(y)
    r    <- y - colMeans(y)                # residuals
    tc   <- rft.critical.t(r, alpha=0.05)  # critical t value
    if ( max(t) > tc ){
        nfp.rft <- nfp.rft + 1
    }
    
    # report results
    msg <- sprintf("Iteration: %d, FPR (RFT): %.3f", n, nfp.rft/n)
    print( msg )
}



