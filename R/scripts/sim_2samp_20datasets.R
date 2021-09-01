
# This script analyzes the 20 datasets in "./Data/20datasets/" using LAAST
#    and saves the critical p-value for each dataset. The results are saved
#    to "./Data/20datasets-laast-results.csv" and are summarized using
#    "./Python/fig_2samp_20datasets.py"
#
# The datasets were generated using the rft1d random number generator
#    (see "./Python/sim_2samp_20datasets.py" for details)


rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


# assemble directories:
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data', '20datasets' )  # path to Data directory
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository


# source the required laast-eval functions:
fnameR1      <- file.path(dirR, 'laast-eval', 'laast.R')
fnameR2      <- file.path(dirR, 'laast-eval', 'loess.R')
fnameR3      <- file.path(dirR, 'laast-eval', 'csv.R')
source( fnameR1 )
source( fnameR2 )
source( fnameR3 )


# process all datasets:
pcrit        <- numeric(20)  # holder for critical p values
false.pos    <- integer(20)  # holder for false positives
for (i in 1:20){
    # load data:
    fname.csv    <- file.path( dirDATA, sprintf( '%d.csv', i ) )
    data.list    <- read.data( fname.csv )
    y1           <- data.list[[1]]
    y2           <- data.list[[2]]
    # run LAAST:
    results      <- mylaast(y1, y2, binSize=2, span=NULL, loess=T, welch=T)
    pc           <- results[[1]]   # critical p value
    df           <- results[[2]]   # data frame
    # check for a false positive:
    fp           <- min(df$pvals) < pc
    false.pos[i] <- fp
    pcrit[i]     <- pc
    print( sprintf('Dataset %d, False positive: %s', i, fp ) )
}
print('-------')
print( sprintf('Total number of false positives: %d', sum(false.pos) ) )


# export results:
fname.out         <- file.path(dirREPO, 'Data', '20datasets-laast-results.csv')
Dataset           <- c(1:20)
Critical.p.values <- format(pcrit, digits=5, scientific=F)
False.positive    <- false.pos
df                <- data.frame( Dataset, Critical.p.values, False.positive )
write.csv(df, file=fname.out, row.names=F, quote=F)

