


rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


# assemble directories 
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data' )  # path to the Data directory in this repository
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository


# source all functions in laast-eval
fileR        <- file.path(dirR, 'LAASTv2', 'LAASTv2.R')
source(fileR)



# load data:
fname.csv    <- file.path(dirDATA, 'Besier2009-MedGastrocF.csv')
data.list    <- read.data( fname.csv )
y1           <- data.list[[1]]
y2           <- data.list[[2]]


# run LAAST:
results      <- mylaast(y1, y2, binSize=2, span=NULL, loess=T, welch=T)
# results      <- mylaast.minimal(y1, y2, alpha=0.05)


# report results:
print( pcrit.laast )  # LAAST-adjusted critical p value (global)
plot( results$q, log(results$pvals), type='l' )

