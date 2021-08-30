
# This script demonstrates how to use the "mylaast" and "mylaast.minimal"
#    functions in ./laast-eval/laast.R
 

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


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
fname.csv    <- file.path(dirDATA, 'Besier2009-MedGastrocF.csv')
# fname.csv    <- file.path(dirDATA, 'Gaussian_FWHM=20_Q=100.csv')
# fname.csv    <- file.path(dirDATA, 'Gaussian_FWHM=20_Q=101.csv')
data.list    <- read.data( fname.csv )
y1           <- data.list[[1]]
y2           <- data.list[[2]]


# run LAAST:
results0     <- mylaast(y1, y2, binSize=2, span=0.17, loess=T, welch=T)
results1     <- mylaast.minimal(y1, y2, alpha=0.05)
pcrit0       <- results0[[1]]
pcrit1       <- results1[[1]]
df0          <- results0[[2]]
df1          <- results1[[2]]


# report results:
print( c(pcrit0, pcrit1) )  # LAAST-adjusted critical p value (global)
plot( df0$q, log(df0$pvals), type='l', col="blue" )
lines( df1$q, log(df1$pvals), col="red" )



