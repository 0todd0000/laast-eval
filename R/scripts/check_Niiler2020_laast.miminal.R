

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


# assemble directories 
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data' )  # path to the Data directory in this repository
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository


# # source all laast-eval functions:
# fnameR1      <- file.path(dirR, 'laast-eval', 'laast.R')
# fnameR2      <- file.path(dirR, 'laast-eval', 'loess.R')
# fnameR3      <- file.path(dirR, 'laast-eval', 'csv.R')
# source( fnameR1 )
# source( fnameR2 )
# source( fnameR3 )


# # load data:
# fname.csv    <- file.path(dirDATA, 'Besier2009-MedGastrocF.csv')
# # fname.csv    <- file.path(dirDATA, 'SimulatedTwoLocalMax.csv')
# # fname.csv    <- file.path(dirDATA, 'Gaussian_FWHM=20.csv')
# data.list    <- read.data( fname.csv )
# y1           <- data.list[[1]]
# y2           <- data.list[[2]]
#
#
# # run LAAST:
# # results      <- mylaast(y1, y2, binSize=2, span=NULL, loess=T, welch=T)
# results      <- mylaast.minimal(y1, y2, alpha=0.05)
#
#
# # report results:
# print( pcrit.laast )  # LAAST-adjusted critical p value (global)
# plot( results$q, log(results$pvals), type='l' )




# source all functions in laast-eval
fileR        <- file.path(dirR, 'LAASTv2', 'LAASTv2.R')
source(fileR)

