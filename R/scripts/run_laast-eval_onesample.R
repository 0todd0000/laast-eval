
# This script demonstrates how to use the "mylaast" and "mylaast.minimal"
#    functions in ./laast-eval/laast.R
 

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


mylaast.onesample <- function(y, binSize=2, alpha=0.05, loess=TRUE, span=NULL) {
    # LAAST implementation
    #
    #     This function can be used to replicate all main LAAST results from Niiler (2020)
    #     "mylaast" is used to distinguish this function from functions in LAASTv2
    #
    #     Comments below detail the connections between this implementation and LAASTv2
    #
    #     The "mlasst.minimal" function above provides a conceptually concise overview of
    #     the proposed LAAST methodlogy, albeit at the expense of substantially reduced
    #     functionality.
    
    # (1) Estimate means, SDs and sample sizes
    if (loess){  # LOESS-smoothed mean and SD estimates
        # estimate span
        #     Optionally skip the computationally-expensive findSpan procedure
        #     by specifying a span value when calling "mylaast"
        if (is.null(span)){
            # LAASTv2 code:  see ./buildFiguresFinal.R Lines 616-618  (in "buildFigure2")
            span  <- findSpan(y)
        }
        # specify x range:
        Q         <- ncol(y)  # numnber of domain nodes
        if (Q%%2==0){  # number of domain nodes is even
            xmin  <- 0
            xmax  <- Q
        } else {
            xmin  <- 1
            xmax  <- Q
        }
        # smooth using LOESS:
        mod   <- mymodelSeries(y, binSize, span, xmin, xmax)
        q     <- mod$times
        m     <- mod$m1
        s     <- mod$sd1
        n     <- mod$n1
        
    } else {  # least-squares estimates 
        q     <- c( 1 : ncol(y) )
        m     <- colMeans( y )
        s     <- apply(y, 2, sd)
        n     <- nrow(y) * rep(1, ncol(y) )
    }

    # (2) Calculate t statistics and conduct inference
    t     <- m / ( s / sqrt(n) )                 # t statistic
    degF  <- n - 1                               # degrees of freedom
    pvals <- 2 * pt( abs(t), degF, lower=FALSE)  # p-values (uncorrected, two-tailed)
    pcrit <- laast.critical.pvalue( t , alpha )  # LAAST's correlation-adjusted critical p value
    df    <- data.frame(q,t,pvals,m,s,n,degF)
    return(   list( pcrit, df)   )
}



# assemble directories:
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data' )  # path to the Data directory in this repository
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository


# source the required laast-eval functions:
fnameR1      <- file.path(dirR, 'laast-eval', 'laast.R')
fnameR2      <- file.path(dirR, 'laast-eval', 'loess.R')
fnameR3      <- file.path(dirR, 'laast-eval', 'csv.R')
source( fnameR1 )
source( fnameR2 )
source( fnameR3 )


# load data:
fname.csv    <- file.path(dirDATA, 'Gaussian_FWHM=20_Q=100.csv')
# fname.csv    <- file.path(dirDATA, 'Gaussian_FWHM=20_Q=101.csv')
data.list    <- read.data( fname.csv )
y            <- data.list[[1]]


# run LAAST:
results      <- mylaast.onesample(y, binSize=2, span=NULL, loess=T)
pcrit        <- results[[1]]
df           <- results[[2]]



# report results:
# matplot( t(y), type='l', lty=1, lwd=0.5, col='blue')
print( pcrit )  # LAAST-adjusted critical p value (global)
plot( df$q, log(df$pvals), type='l', col="blue", ylim=c(-5,1) )
abline(h=log( pcrit ), col='red', lty=3)


