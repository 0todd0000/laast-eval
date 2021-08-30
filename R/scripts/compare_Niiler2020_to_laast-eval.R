

# This script checks whether the results from "mylaast" in ./laast-eval/laast.R
#    yields the same results as LAASTv2, using the Fig.2 data from Niiler (2020)
#
# NOTE: 
#    There appears to be an error in LAASTv2 on Line 1189 of buildFiguresFinal.R:
#       The call to "smoothAndFit" uses "type", but I believe that it should use "type1"
#
#       This causes smoothAndFit to uses the T value for the Welch test adjustment:
#           (Line 624):  t = calcT(m1,m2,se1a,se2a,n1a,n2a);
#   
#       BUT to use the P value for the Pooled test adjustment:
#           (Line 661):  t = calcPTpooled(m1,m2,sp,n1a,n2a);
#
#       I think this is simply a coding error, as a I see no rationale for
#       using T-values for one case and P-values for the other.
#
#       Regardless, the current version of laast-eval was written to replicate
#       the results of LAASTv2, so like LAASTv2 it uses T-values in one case
#       and P-values in the other.



rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


# assemble directories 
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data' )  # path to the Data directory in this repository
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository


#source all laast-eval functions:
fnameR1      <- file.path(dirR, 'laast-eval', 'laast.R')
fnameR2      <- file.path(dirR, 'laast-eval', 'loess.R')
fnameR3      <- file.path(dirR, 'laast-eval', 'csv.R')
source( fnameR1 )
source( fnameR2 )
source( fnameR3 )


# load results from Niiler (2020)
#    these results were calculated and exported using "Niiler2020_fig2_export.R"
fnamePCRIT  <- file.path( dirDATA, 'Niiler2020-results', 'fig2-pcrit.csv')
fname0      <- file.path( dirDATA, 'Niiler2020-results', 'fig2-pooled.csv')
fname1      <- file.path( dirDATA, 'Niiler2020-results', 'fig2-welch.csv')
dfPCRIT     <- read.csv( fnamePCRIT )
df0         <- read.csv( fname0 )
df1         <- read.csv( fname1 )
pcrit0      <- dfPCRIT[1,1]
pcrit1      <- dfPCRIT[1,2]


# calculate results using laast-eval
fname.csv    <- file.path(dirDATA, 'Besier2009-MedGastrocF.csv')
data.list    <- read.data( fname.csv )
y1           <- data.list[[1]]
y2           <- data.list[[2]]
# resultsMIN   <- mylaast.minimal(y1, y2, alpha=0.05)
span         <- NULL   # setting span to NULL will trigger mylaast to estimate it
span         <- 0.17   # optionally use this to bypass the computationally expensive span estimate
results0     <- mylaast(y1, y2, binSize=2, span=span, loess=T, welch=F)  # pooled
results1     <- mylaast(y1, y2, binSize=1, span=span, loess=T, welch=T)  # Welch (./LAASTv2/buildDiguresFinal.R uses binSize=1 on Line 663)
pcritEVAL0   <- results0[[1]]
pcritEVAL1   <- results1[[1]]
dfEVAL0      <- results0[[2]]
dfEVAL1      <- results1[[2]]
# print( c(pcritEVAL0,pcritEVAL1) )


# plot:
graphics.off()


# pdf('/Users/todd/Desktop/fig.pdf')
plot( df0$times, df0$logp, type='l', col='blue', xlab="Time (%)", ylab="log(p)", ylim=c(-15,0) )
lines( df1$times, df1$logp, col='red' )
lines( dfEVAL0$q, log( dfEVAL0$pvals ), col='blue', lty=3, lwd=5 )
lines( dfEVAL1$q, log( dfEVAL1$pvals ), col='red', lty=3, lwd=5 )
col0 <- rgb(0.2, 0.8, 1)
col1 <- rgb(1, 0.5, 0.2)
abline(h=log( pcrit0 ), col=col0)
abline(h=log( pcritEVAL0 ), col=col0, lty=3, lwd=5)
abline(h=log( pcrit1 ), col=col1)
abline(h=log( pcritEVAL1 ), col=col1, lty=3, lwd=5)
labels <- c("Niiler (2020) - pooled", "Niiler (2020) - Welch", "laast-eval - pooled", "laast-eval - Welch")
labels <- c(labels, "", "[pcrit] Niiler (2020) - pooled", "[pcrit] Niiler (2020) - Welch", "[pcrit] laast-eval - pooled", "[pcrit] laast-eval - Welch")
legend(30, -10, legend=labels, col=c("red", "blue", "red", "blue", "white", col0, col1, col0, col1), lty=c(1,1,3,3,0,1,1,3,3), lwd=c(1,1,5,5,0,1,1,5,5), cex=0.8)
title()

