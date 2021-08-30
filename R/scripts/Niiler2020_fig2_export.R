

rm( list = ls() ) # clear workspace
graphics.off()    # close all graphics


buildFigure2_tp <- function(inputData,RFTtvalues) {
    # {TP} This is a modification of the "buildFigure2" from "./LAASTv2/buildFiguresFinal.R"
    #      The purpose of the modifications is to faciliate results exporting.
    #      The "p1" and "p2" variables are now returned so that relevant data can be exported.
    #      Otherwise this function remains essentially unaltered.
    #
    #      All source code in ./LAASTv2/ remains unaltered.
    #
    #      Modifications made by Todd Pataky on 2021-08-30
    #      See the original "buildFigure2" function for more details.

    RFTinfo = c(2.777,41);
    N1 = 27;	# Number in each group
    N2 = 16;

    # Find optimal smoothing parameter
    ind1 = which(inputData$Group == 'Control');
    ind2 = which(inputData$Group == 'Pain');
    x1 = inputData$Time[ind1];
    x2 = inputData$Time[ind2];
    y1 = inputData$MedGastrocF[ind1];
    y2 = inputData$MedGastrocF[ind2];

    #

    s1 = findSpan(x1,y1);
    s2 = findSpan(x2,y2);
    span = min(c(s1,s2));
    print(paste('span = ',span));


    # Set up graphical parameters
    tRa = c(1,6.5);	# Ranges of t, dof, and p plots
    dfRa = c(20,42);
    pRa = c(-18,0);

    showX = c(F,F,F,T);	# which plots to show scales and legends on
    showY = c(T,T,T,T);
    showL = c(F,F,F,T);
    lPosMean = c(0.2,0.7);	# legend positions
    #lPosP = c(0,0);		
    lPosP = c(0.5,0.13);
    #lPosP = c(0.1,0.35);	
    showPts = F;		# show points on mean plot or no
    lHoriz = T;	# horizontal p legend

    # Graphical parameters for LAAST
    gParams = list(0,100,20,"% Cycle","Musle Force (N)",span,"MedGast",tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz)

    # Do pooled t-test
    p1 = laastVsRFT(inputData,'Time','MedGastrocF',c('Pain','Control'),2,'pp',gParams,RFTinfo,RFTtvalues);
    alphaN1 <- alphaN

    # p1a = p1[[1]]+annotate("text",x=-3,y=1150,label="(a)");
    # p2a = p1[[2]]+annotate("text",x=-3,y=40,label="(c)");
    # p3a = p1[[3]]+annotate("text",x=-3,y=5.75,label="(b)")+
    # annotate("text",x=20,y=2.5,label="LAAST",color="#10BE44")+
    # annotate("text",x=10,y=4.5,label="RFT",color="blue");
    # p4a = p1[[4]]+annotate("text",x=-3,y=-1,label="(d)");

    showX = c(F,F,F,T);
    showY = c(F,F,F,F);
    showL = c(T,F,F,T);
    lPosMean = c(0,0);
    lPosP = c(0.5,0.13);
    showPts = F;
    lHoriz = T;	# horizontal p legend

    # Graphical parameters for LAAST
    gParams = list(0,100,20,"% Cycle","Knee Flexion (deg)",span,"MedGast",tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz)

    # Do Welch t-test
    p2 = laastVsRFT(inputData,'Time','MedGastrocF',c('Pain','Control'),1,'p',gParams,RFTinfo,RFTtvalues);
    alphaN2 <- alphaN

    # p1b = p2[[1]]+annotate("text",x=-3,y=1150,label="");
    # p2b = p2[[2]]+annotate("text",x=-3,y=40,label="(f)");
    # p3b = p2[[3]]+annotate("text",x=-3,y=5.75,label="(e)")+
    #         annotate("text",x=10,y=5.5,label="LAAST",color="#10BE44")+
    #         annotate("text",x=20,y=2.5,label="RFT",color="blue");;
    # p4b = p2[[4]]+annotate("text",x=-3,y=-1,label="(g)");
    #
    # # Make a blank plot
    # df1 = data.frame();
    # pb = ggplot(df1) + geom_point() +
    #     xlim(0,5)+ ylim(0,5) +
    #     theme(line=element_blank(),text=element_blank(),panel.background = element_rect(fill = "white"))

    # g1 <- ggarrange(p1a,pb, p3a,p3b, p2a,p2b, p4a,p4b, nrow =4);    # from egg - better for
    # print(g1);

    # fN = "Figure2.png";
    # ggsave(filename=fN,plot=g1,width=8,height=8,dpi=400);  # units="in",
    return( list(alphaN1, alphaN2, p1, p2) )
}




# assemble directories 
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data' )  # path to the Data directory in this repository
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository
dirLAASTv2   <- file.path(dirR, 'LAASTv2')


# source all needed R files:
fileR1       <- file.path(dirLAASTv2, 'LAASTv2.R')
fileR2       <- file.path(dirLAASTv2, 'buildFiguresFinal.R')
source(fileR1)
source(fileR2)
setwd( dirLAASTv2 )  # set the working directory so R knows where to find the LAASTv2 data files


# run LAAST  (takes approximatel 5-10 s to run and create all figures)
inputData  <- readPatakyData()
RFTtvalues <- readPatakyT()
results    <- buildFigure2_tp( inputData, RFTtvalues )


# assemble and export critical p values:
pcrit_pooled <- results[[1]]  # critical p value (pooled SD)
pcrit_welch  <- results[[2]]  # critical p value (Welch)
df           <- data.frame( pcrit_pooled, pcrit_welch )
fnameCSV     <- file.path( dirDATA, 'Niiler2020-results', 'critical_p_values.csv')
write.csv( df, fnameCSV, row.names=F)


# assemble and export pooled results ( t and p values ):
p1       <- results[[3]]
times    <- p1[[3]][[1]]$times
t        <- p1[[3]][[1]]$sig
logp     <- p1[[4]][[1]]$sig
df       <- data.frame(times, t, logp)
fnameCSV <- file.path( dirDATA, 'Niiler2020-results', 'fig2-pooled.csv')
write.csv( df, fnameCSV, row.names=F)


# assemble and export Welch results ( t and p values ):
p2       <- results[[4]]
times    <- p2[[3]][[1]]$times
t        <- p2[[3]][[1]]$sig
logp     <- p2[[4]][[1]]$sig
df       <- data.frame(times, t, logp)
fnameCSV <- file.path( dirDATA, 'Niiler2020-results', 'fig2-welch.csv')
write.csv( df, fnameCSV, row.names=F)




