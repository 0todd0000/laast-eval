

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


read.Besier.data <- function(file.name){
    a         <- read.csv(fname.csv, header=FALSE)
    a         <- as.matrix( a )
    group     <- a[,1]
    y         <- a[,2:ncol(a)]
    u         <- unique(group)
    y1        <- y[group==u[1] , ]   # first group
    y2        <- y[group==u[2] , ]   # second group
    return( list(y1, y2) )
}

# load data:
fname.csv    <- 'Besier2009-MedGastrocF.csv'
data.list    <- read.Besier.data( fname.csv )
y1           <- data.list[[1]]
y2           <- data.list[[2]]

# source all files in laast-eval
dirR         <- dirname( dirname( sys.frame(1)$ofile ) )  # path to the R directory in this repository
dirLAASTeval <- file.path(dirR, 'laast-eval')
filesR       <- list.files(dirLAASTeval, pattern="*.R", full.names=TRUE)
sapply(filesR, source)



binSize = 2
# results = mylaast(y1, y2, binSize)
results = mylaast(y1, y2, binSize, span=NULL, loess=T, welch=T)
# results = mylaast.minimal(y1, y2, alpha=0.05)



print( pcrit.laast )  # LAAST-adjustuded critical p value (global)
# plot( results$pvals )
plot( results$q, log(results$pvals), type='l' )



