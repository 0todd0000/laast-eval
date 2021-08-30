
laast.critical.pvalue <- function(t, alpha=0.05){
    # Correlation-based critical p-value adjustment
    #
    # This is the core of the proposed LAAST method.
    #     This adjustment does not appear to be theoretically
    #     validated, and numerical simulations show that it
    #     can not control false positive rates.
    #     USE OF THIS ADJUTMENT IS NOT RECOMMENDED 
    #
    # Variable notes:
    #     In LAASTv2, and for the Welch case, this adjustment
    #     uses "calcPTpooled", which returns P-values and not
    #     T-values; thus "t" below actually refers to P-values
    #     for the Welch case.
    #
    #     Regardless, since the P and T value functions exhibit
    #     approximately the same moothness, using P vs. T is not
    #     expected to markedly change the correlation strength "rho"
    t      <- t[!is.na(t)]
    tL     <- t[1:length(t)-1]
    tR     <- t[2:length(t)]
    rho    <- cor(tL, tR)
    alphaN <- alpha / ( length(t)*(1-rho^2)+ rho^2 )
    return( alphaN )
}





mylaast.minimal <- function(y1, y2, alpha=0.05) {
    # Minimal LAAST implementation
    #
    #    This implementation uses least-squares estimates of
    #    means and standard deviations, and also a simple
    #    T value calculation that is based on the pooled
    #    standard deviation.
    #
    #    This function was prepared to emphasize that LAAST's
    #    main proposed novelty is the alpha-adjustment as
    #    implemented in laast.critical.pvalue
    #
    #    This function was also prepared for numerical
    #    validation purposes; it is considerably faster
    #    than using LOESS smoothing. Moreo
    #
    # means, SDs and sample sizes (least-squares estimates)
    q     <- c( 1 : ncol(y1) )
    m1    <- colMeans( y1 )
    m2    <- colMeans( y2 )
    s1    <- apply(y1, 2, sd)
    s2    <- apply(y2, 2, sd)
    n1    <- nrow(y1) * rep(1, ncol(y1) )
    n2    <- nrow(y2) * rep(1, ncol(y2) )
    # calculate t statistic and degrees of freedom (using pooled SD)
    sp    <- sqrt(     ( (n1-1)*s1*s1 + (n2-1)*s2*s2 ) / (n1 + n2 - 2)     )
    t     <- (m1 - m2) / (sp * sqrt( 1/n1 + 1/n2 ) )
    degF  <- n1 + n2 - 2  # degrees of freedom
	# calculate (uncorrected, two-tailed) p values
	pvals <- 2 * pt( abs(t), degF, lower=FALSE)
    # calculate the LAAST correlation-adjusted critical p value
    pcrit <- laast.critical.pvalue( pvals, alpha )
    df    <- data.frame(q,t,pvals,m1,m2,s1,s2,n1,n2,degF)
    return(   list( pcrit, df)   )
}



mylaast <- function(y1, y2, binSize=2, alpha=0.05, loess=TRUE, span=NULL, welch=TRUE) {
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
            span1 <- findSpan(y1)
            span2 <- findSpan(y2)
            span  <- min( c(span1,span2) )
        }
        # specify x range:
        # LAASTv2 code: see ./buildFiguresFinal.R Lines 660 (in "buildFigure2")
        #     xmin and xmax are specified manually there as 0 and 100, respectively
        #     Here xmin and xmax are calculated algorithmically.
        Q         <- ncol(y1)  # numnber of domain nodes
        if (Q%%2==0){  # number of domain nodes is even
            xmin  <- 0
            xmax  <- Q
        } else {
            xmin  <- 1
            xmax  <- Q
        }
        # smooth the first group using LOESS:
        # LAASTv2 code:  see ./LAASTv2.R Lines 590 and 597 (in function "smoothAndFit")
        #   This LOESS smoothing procedure is implemented in ./laast-eval/loess.R
        mod1  <- mymodelSeries(y1, binSize, span, xmin, xmax)
        q     <- mod1$times
        m1    <- mod1$m1
        s1    <- mod1$sd1
        n1    <- mod1$n1
        # smooth the second group using LOESS:
        mod2  <- mymodelSeries(y2, binSize, span, xmin, xmax)
        m2    <- mod2$m1
        s2    <- mod2$sd1
        n2    <- mod2$n1
        
    } else {  # least-squares estimates 
        # LAASTv2 code:  this is not implemented in LAASTv2
        #     It is implemented here just to show that the proposed inference procedure
        #     can also be used with non-smoothed data.
        q     <- c( 1 : ncol(y1) )
        m1    <- colMeans( y1 )
        m2    <- colMeans( y2 )
        s1    <- apply(y1, 2, sd)
        s2    <- apply(y2, 2, sd)
        n1    <- nrow(y1) * rep(1, ncol(y1) )
        n2    <- nrow(y2) * rep(1, ncol(y2) )
    }

    # (2) Calculate t statistic and degrees of freedom
    if (welch){  # Welch test with # Welch-Satterthwaite degrees-of-freedom correction
        # LAASTv2 code:  see ./LAASTv2.R Lines 1153 and 1125 (in functions "calcT" and "calcDF")
        t     <- (m1 - m2) / sqrt(   (s1^2)/n1 + (s2^2)/n2   )
        v1    <- n1-1
        v2    <- n2-1
        degF  <- ( (s1^2)/n1 + (s2^2)/n2)^2/( (s1^4)/(n1*n1*v1) + (s2^4)/(n2*n2*v2))
        
    } else {  # statistic using pooled SD
        # LAASTv2 code:  see ./LAASTv2.R Lines 939 and 1182 (in functions "getPooledSD" and "calcTpooled")
        sp    <- sqrt(     ( (n1-1)*s1*s1 + (n2-1)*s2*s2 ) / (n1 + n2 - 2)     )
        t     <- (m1 - m2) / (sp * sqrt( 1/n1 + 1/n2 ) )
        degF  <- n1 + n2 - 2  # degrees of freedom
    }
    
    # (3) Calculate (uncorrected, two-tailed) p values
    # LAASTv2 code:  see ./LAASTv2.R Lines 1141 and 1174 (in functions "calcPT" and "calcPTpooled")
    pvals <- 2 * pt( abs(t), degF, lower=FALSE)
    
    # (4) Calculate LAAST's correlation-adjusted critical p value
    if (welch){
        # LAASTv2 code:  see ./LAASTv2.R Lines 624-637 (in function "smoothAndFit")
        pcrit <- laast.critical.pvalue( t, alpha )
    } else {
        # LAASTv2 code:  see ./LAASTv2.R Lines 661-674 (in function "smoothAndFit")
        #    Unlike the adjustment for the Welch case above, which uses T-values,
        #    this adjustement uses P-values (note that "calcPTpooled" returns P-values)
        #    This is likely a coding mistake in LAASTv2, but is replicated here so that
        #    identical results can be obtained from laast-eval and LAASTv2.
        pcrit <- laast.critical.pvalue( pvals, alpha )  # Bizarre
    }
    
    df    <- data.frame(q,t,pvals,m1,m2,s1,s2,n1,n2,degF)
    return(   list( pcrit, df)   )
}





