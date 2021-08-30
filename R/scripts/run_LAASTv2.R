
# This script uses LAASTv2 to reproduce Fig.2 from Niiler (2020)


rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


# assemble directories 
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirLAASTv2   <- file.path(dirREPO, 'R', 'LAASTv2')


# source all needed R files:
fileR1       <- file.path(dirLAASTv2, 'LAASTv2.R')
fileR2       <- file.path(dirLAASTv2, 'buildFiguresFinal.R')
source(fileR1)
source(fileR2)
setwd( dirLAASTv2 )  # set the working directory so R knows where to find the LAASTv2 data files


# run LAAST  (takes approximatel 5-10 s to run and create all figures)
inputData  <- readPatakyData()
RFTtvalues <- readPatakyT()
results    <- buildFigure2( inputData, RFTtvalues )
