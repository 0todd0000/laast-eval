



findSpan <- function(y) {
    
    # Store these as globals to evaluate smoothing
    # {TP}:  these three variables are defined outside this function in LAASTv2,
    #        but it appears that they can be confined to the local namespace
    #        without affecting the results.
    rmsArr = c(); # RMS values found in the smoothing iterations
    rmsCnt = 1;	  # index of rmsArr and spanArr
    spanArr = c();# span values found in the smoothing iterations
    
    
	# WARNING: This takes a while to run.
	#
	# Find smoothing span by k-fold cross-validation that finds minimum RMS 
	# on leave-one out data.  Please see:
	# J.S. Lee, D.D. Cox, Robust Smoothing: Smoothing Parameter Selection 
	#	and Applications to Fluorescence Spectroscopy, Comput Stat Data Anal. 54 
	#	(2010) 3131–3143. doi:10.1016/j.csda.2009.08.001.

	# x is a 1d vector of x values and 
	# y is the corresponding 1d vector of y values

	# Outputs s0, the smoothing parameter or span
	
    # create X array
    x  <- c( 1 : ncol(y) )
    X  <- t( replicate( nrow(y), x ) )
    x  <- as.vector( t(X) )
    y  <- as.vector( t(y) )
    
    
	print('Looking for best smooth... this may take a while...');
	dfXY = data.frame(x,y);
	
	# Sort the data prior to all this
	dfXY1  <- dfXY[order(dfXY$x),] 	
	x = dfXY1$x;
	y = dfXY1$y;
	
	# In loess, spans define the amount of the data that is used for smoothing.
	# A span of 1 uses all the data, while a span of 0.1 uses windowing on 1/10
	# of the data at a time.  This is a parameter that should be adjusted dependent
	# on the specific data.  Spectroscopy will use smaller spans, in general, since
	# emission or absorption lines will be relatively sharp deviations from the
	# rest of the curve.  In biomechanics, most spans will be much larger since,
	# except for impact forces, the curves will be smoother.
	
	span0 = 0.05;	# Starting span 
	dSpan = 1.0 - span0;	# difference from max
	dS = 0.02;		# How much to adjust the span upward each iteration
	N = floor(dSpan/dS); # Number of iterations
	
	#N = 39;
	maxCnt = 5;	# Speed things up if no better value is found for a while
	
	s0 = span0;
	rmsLast = 1e8;
	
	Nsplits = 5;	# Number of folds or splits in k-fold cross-validation scheme

	rmsArr <<- c();	# Reinitialize
	spanArr <<- c();
	rmsCnt <<- 1;

	for (j in c(1:Nsplits) ) {
	
		cnt = 0;	

		nx = length(x);
	
		a = c(1:nx);		
		# sample every Nsplits points from a
		aS <- a[seq(j, length(a), Nsplits)]

		xS <- x[-aS];	# Leave points with aS indices out
		yS <- y[-aS];	
		xS2 <- x[aS];	# Points that are left out
		yS2 <- y[aS];
		
		for (i in c(1:N) ) {
			
			#print(paste('Current span = ',span0));
			sp0 = loess(yS~xS,span=span0);
			
			sp0P = predict(sp0,xS2,se=TRUE);
			m1 = sp0P$fit;	# mean of model to data at values xS2

			res = m1 - yS2;	# Residuals

			n = length(res);
		
			res2 = res*res
			rms = sqrt( sum( res2 ) / n );
			rmsArr[rmsCnt] <<- rms;
			spanArr[rmsCnt] <<- span0;
			rmsCnt <<- rmsCnt + 1; 
			#print(paste("RMS = ",rms));
			if (rms < rmsLast) {
				rmsLast = rms;
				s0 = span0;
				print(paste("Current best span is = ",s0));
				span0 = span0+dS # add extra dS to speed things up.
				cnt = 0;
			} 
			
			span0 = span0 + dS;	# Update the smoothing parameter (again, as necessary)
			if (cnt == maxCnt) break;
			if (span0 > 1) break;
			
			cnt = cnt+1;
		}
	
	}
	print("OK, smoothing parameter found...");
	return (s0);
	
}





nPerSamp <- function(timePts,N,xmin,xmax) {
	# Estimate the number of points at each sample along the time domain.
	# This function uses a density function to model the frequencies
	# of data across all time points.  This density estimate is then
	# resampled at N intervals and scaled up to the total number of
	# original data points.  
	
	# This is a sampling approach to the data rather than a binning approach,
	# where the number of samples smoothly changes throughout the time
	# domain, and so, any sample will have the same approximate properties
	# of any sample around it.
	
	# Extend timePts to avoid edge effects
	#x1 = timePts-max(timePts);
	#x2 = timePts+max(timePts);
	#x = c(x1,timePts,x2);
	timePts = sort(timePts);	# has to be in order or this doesn't work	
	
	b = -timePts;
	dba = min(timePts)-max(b);
	rTPleft = -rev(timePts)+dba-1
	
	dbb = max(timePts)-min(b);
	rTPright = rev(-timePts)+dbb+1
	
	x = c(rTPleft,timePts,rTPright);
	#
	
	d = density(x,na.rm=TRUE);
	
	d1x = d$x;
	d1y = d$y
	n = 3*length(timePts);
	d1y = n*d1y; # scale to total number of points
	#plot(d);
	
	ind = which(d1x > xmin & d1x < xmax);
	#ind = which(d1x > 0 & d1x <101);
	d1x = d1x[ind];
	d1y = d1y[ind];
	 
	appr = approx(d1x,d1y,n=N);  # Resampling or approximation
	#appr = spline(d1x,d1y,N); # Now there's an endpoint problem
	d2y = appr$y;	# undercounts a touch? - Feb 1, 2018
	#print(paste("N = ",mean(d2y)))
	#appr = sampleEvery(d1y,N);
	#print(appr);
	#print(paste("N samp = ",N));
	#print(appr)
	#d1y = appr
	 
	return (d2y);
	
}


mymodelSeries <- function(y, binSize0, span0, xmin=0, xmax=100) {
	
	# binSize0 is the bin size along the x axis.  This determines how much data
	# 	are grouped together for analysis.
	
	library(msir)

    # create X array
    x  <- c( 1 : ncol(y) )
    X  <- t( replicate( nrow(y), x ) )
    X  <- as.vector( t(X) )
    Y  <- as.vector( t(y) )

    # xmin  = min(x)
    # xmax  = max(x)
    nsamp = (xmax-xmin)/binSize0;
    times = c(1:nsamp)*binSize0+ xmin;	# Well behaved x values for model

    # df0   = data.frame(x,y);
    # sp0   = loess(y~x,df0,span=span0)  # smoothed using loess
    df0   = data.frame(X, Y)
    sp0   = loess(Y~X, df0, span=span0)  # smoothed using loess
    sp0P  = predict(sp0, times, se=TRUE);  # central (i.e., mean) estimate
    m1    = sp0P$fit;	# mean of model to data at values xArr

    # Get standard deviation of the data
    # lsd = loess.sd(x,y,nsigma=1,span=span0);
    lsd = loess.sd(X, Y, nsigma=1,span=span0);
    y1 = lsd$sd;
    x1 = lsd$x;
    df1 = data.frame(x1,y1);
    sp1 = loess(y1~x1,df1,span=span0);
    sp1P = predict(sp1,times,se=TRUE);
    sd1 = sp1P$fit;  # modeled standard deviation of data at values xArr

    # Estimate number of data points at each sample based on density of x
    n1 = nPerSamp(X,nsamp,xmin,xmax);
    # #n1 = nPerBin(x,1,100,binSize0);

    model <- data.frame(times,m1,n1,sd1)
    return (model)	
}



