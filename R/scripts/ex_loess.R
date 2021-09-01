
# This script demonstrates a potential biomechanical problem
#    with the LOESS procedure as implemented in LAASTv2.
#
#    The data are ground reaction forces (GRF) during fast
#    fast walking. As seen in the results, the LOESS procedure
#    in LAASTv2 applies a constant smoothing parameter across
#    the whole time domain, which disrupts the high frequency
#    GRF content near impact (i.e., shortly after time=0).
#    This may be non-ideal biomechanically.


rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


# assemble directories:
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data' )  # path to the Data directory in this repository
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository


# source the required laast-eval functions:
fnameR       <- file.path(dirR, 'laast-eval', 'loess.R')
source( fnameR )


# load data:
fnameCSV     <- file.path( dirDATA, 'Pataky2008-grf.csv')
y            <- as.matrix( read.csv(fnameCSV, header=FALSE) )


# run LOESS
# span        <- findSpan(y)
span        <- 0.19   # optionally set "span" value to avoid repeating the "findSpan" calculation
xmin        <- 0
xmax        <- ncol(y)-1
binSize     <- 2
mod1        <- mymodelSeries(y, binSize, span, xmin, xmax)
q           <- mod1$times
m           <- mod1$m1


# plot results:
matplot( t(y), type='l', lty=1, lwd=0.5, col='blue')
lines(q, m, col='red', lwd=3)


# export results:  (paper figure created using Python; see: "fig_loess.py")
fnameOUT    <- file.path( dirDATA, 'Pataky2008-grf-loess-mean.csv')
a           <- rbind( q, m )
write.table(format(a, digits=3, scientific=F), file=fnameOUT, sep=',', row.names=F, col.names=F, quote=F)




