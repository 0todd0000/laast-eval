
# This script calculates the false positive rates in one-sample tests for:
#    - LAAST  (see "mylaast.onesample" in ./R/laast.R)
#    - RFT (random field theory;  see "rft.critical.t" in ./R/sim.R)
 

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics



# (0) Source R files:
dirREPO  <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirR     <- file.path( dirREPO, 'R' )     # path to the R directory in this repository
# source the required laast-eval functions:
fnameR1  <- file.path(dirR, 'laast-eval', 'laast.R')
fnameR2  <- file.path(dirR, 'laast-eval', 'loess.R')
fnameR3  <- file.path(dirR, 'laast-eval', 'sim.R')
source( fnameR1 )
source( fnameR2 )
source( fnameR3 )




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
nfp.laast <- 0   # number of false positives (LAAST)

# run simulation:
set.seed(2)
for (i in 1:1000){
    n    <- n + 1                          # increment total number of datasets
    y    <- rnorm1d(J, Q, nbasis, norder)  # generate random dataset
    
    # RFT inference:
    t    <- onesample.t.stat(y)            # typical, least-squares, one-sample t statistic
    r    <- y - colMeans(y)                # residuals
    tc   <- rft.critical.t(r, alpha=0.05)  # critical t value (calculated from the residuals of y)
    if ( max(t) > tc ){                    # false positive encountered
        nfp.rft <- nfp.rft + 1             # increment number of false positives
    }
    
    # LAAST inference:
    results  <- mylaast.onesample(y, binSize=2, span=NULL, loess=T)
    pc       <- results[[1]]         # critical p value
    p        <- results[[2]]$pvals   # p values
    if ( min(p) < pc ){              # false positive encountered
        nfp.laast <- nfp.laast + 1   # increment number of false positives
    }

    # report results
    msg <- sprintf("Iteration: %d, FPR (RFT): %.3f, FPR (LAAST): %.3f", n, nfp.rft/n, nfp.laast/n)
    print( msg )
}



