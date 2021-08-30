# buildDistributions.R

## This code is what was used to produce the paper:
## 	Comparing Groups of Time Dependent Data Using LOESS Alpha-Adjusted Serial T-tests ##
## 		by Timothy Niiler, 
##  	Department of Physics, Penn State Brandywine
##  	Gait Laboratory, Nemours/AI duPont Hospital for Children
##

#	Copyright Timothy Niiler 2015-2019
#
#      This program is free software; you can redistribute it and/or modify
#      it under the terms of the GNU General Public License as published by
#      the Free Software Foundation; either version 2 of the License, or
#      (at your option) any later version.
#      
#      This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#      
#      You should have received a copy of the GNU General Public License
#      along with this program; if not, write to the Free Software
#      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#      MA 02110-1301, USA.

###### Creation of Tables C.1 and C.2 in Appendix C.

# C.1 
#	Generate the distribution using dist1 = buildNullDistribGaussianNoise(10000)
#		This distribution may also be generated using dist1 = buildDistribNonNorm(10000,alpha), where alpha = 0
#
# 	Then use sampleNullTest() as shown below to test LAAST
#		out5 = sampleNullTest(N1=5,N2=5, d=0, dfIn=dist1, Nrep=1000);
#		d is the effect size
#		Nrep is the number of repetitions
#			If N1=N2 is larger, reduce this unless your computer is very fast or you have lots of time.
#			For N1=N2=5, Nrep=1000.
#		N1 is the sample size for the first group
#		N2 is the sample size for the second group
#		dfIn is the data frame generated previously.
#	rft1d may be run via the Python script: nullTestRFT.py.

# C.2
#	Generate the distribution using dist2 = buildDistribNonNorm(10000,alpha), where alpha = 0, 1, or 2
#
# 	Then use sampleNullTest() as shown below to test LAAST
#		out5 = sampleNullTest(N1=5,N2=5, d=0, dfIn=dist2, Nrep=1000);
#		d is the effect size (0.5, 1.0, 1.5)
#		Nrep is the number of repetitions
#			If N1=N2 is larger, reduce this unless your computer is very fast or you have lots of time.
#			For N1=N2=5, Nrep=1000.
#		N1 is the sample size for the first group
#		N2 is the sample size for the second group
#		dfIn is the data frame generated previously.
#	rft1d may be run via the Python script: nullTestRFTskew.py.


# Last edit: May 23, 2019
buildNullDistribGaussianNoise <- function(Ncurves) {
	# Used for Table C.1
	# Add Gaussian noise to the signal as an offset (mean of zero, sd of 1);
	
	X = c(1:360);
	Xrad = (pi/180)*X;
	Y = sin(Xrad)+rnorm(1,0,1);
	Subject = rep(0,360);
	
	
	for (i in c(1:Ncurves) ) {
		X1 = c(1:360);
		X1rad = (pi/180)*X1;
		Y1 = sin(X1rad)+rnorm(1,0,1);
		Subject1 = rep(i,360);
	
		
		X = c(X,X1);
		Y = c(Y,Y1);
		Subject = c(Subject,Subject1);
		
	}

	dfOut = data.frame(X,Y,Subject);
	
	# Visualize
	library(ggplot2);
	
	p1 = ggplot(dfOut,aes(x=X,y=Y))+geom_point(alpha=0.2);
	
	print(p1);
	
	return (dfOut);
}

# Last edit: May 28, 2019
buildDistribNonNorm <- function(Ncurves,alpha) {
	# Used for Table C.2	
	# Add Gaussian noise to the signal as an offset (mean of zero, sd of 1);
	
	X = c(1:360);
	Xrad = (pi/180)*X;
	Y = sin(Xrad)+randn_skew_fast(1,alpha,loc=0.0,scale=1.0);
	
	Subject = rep(0,360);	
	
	for (i in c(1:Ncurves) ) {
		X1 = c(1:360);
		X1rad = (pi/180)*X1;
		#Y1 = sin(X1rad)+rnorm(1,0,1);
		Y1 = sin(X1rad)+randn_skew_fast(1,alpha,loc=0.0,scale=1.0);
		Subject1 = rep(i,360);
	
		
		X = c(X,X1);
		Y = c(Y,Y1);
		Subject = c(Subject,Subject1);
		
	}

	dfOut = data.frame(X,Y,Subject);
	
	# Visualize
	library(ggplot2);
	
	p1 = ggplot(dfOut,aes(x=X,y=Y))+geom_point(alpha=0.2);
	
	print(p1);
	
	return (dfOut);
}


rateArr <<- c();
rcnt <<- 1;

# Last edit: May 28, 2019
# Instead of running sampleNullTest a bunch of separate times, batch it!
# Warning, this will take a long time, so go get dinner, catch a movie, or two,
# take a nap, etc.  
batchNullTest <- function(dfIn,alpha) {
	# dfIn = input data created using buildDistribNonNorm()
	#	has form:  data.frame(X,Y,Subject)
	# alpha = the skew factor (aka gamma in paper).  This is only used to keep track of 
	#	alpha for the output file which is used by Python.
	#
	# Note that this function is used only to get an estimate of the effects.  Not enough 
	#	iterations are done for real accuracy in estimation.
	#
	outString = '';
	dVec = c(0,0.5,1,1.5);	# Effect sizes to test
	nTimes = 100;
	fn = paste('batchOutputAlpha',alpha,'.txt',sep='');
	sink(fn, append=TRUE, split=TRUE);
	for (i in c(1:3) ) {
		d = dVec[i];
		# These have the form sampleNullTest(N group 1, N group 2, effect size, dataIn, N iterations or draws)
		out5 = sampleNullTest(5,5, d, dfIn, 1000);
		out10 = sampleNullTest(10,10, d, dfIn, 250);
		out20 = sampleNullTest(20,20, d, dfIn, 100);
		
		outString = paste(outString,out5,out10,out20,':::');

	}
	sink();
	
	print(outString);
	return (outString);
}


# Last edit: May 28, 2019
sampleNullTest <- function(N1,N2,d, dfIn, Nrep) {
	# Find the false discovery rate for LAAST given 
	# a null distribution
	# N1 = number of samples drawn from distribution
	# N2 = ""
	# d = effect size = (mu2-mu1)/sd, use zero for null distrib
	# dfIn = data frame containing null distribution
	# Nrep = how many times to sample
	
	subject = dfIn$Subject;
	us = unique(subject);
	
	pvalArr = c();
	fp = c();
	rateArr <<- c();	# Reinitialize
	rcnt <<- 1;
	
	if (d != 0) {	# Otherwise, it's wasted memory
		dfCopy = dfIn
		dfCopy$Y = dfCopy$Y + d;	# copy with offset for power test
	}
	
	for (i in c(1:Nrep)) {
		print(paste("**************** DRAW Number ",i," ****************"));
		# Sampling from distribution (specifically, pick N1 and N2 subjects with replacement)
		s1 = sample(us,N1,replace=TRUE);
		s2 = sample(us,N2,replace=TRUE);

		# Get the data for first group based on sampled subject matches
		ind1 = c();
		for (s in s1) {
			ind = which(dfIn$Subject == s);
			ind1 = c(ind1,ind);
		}

		# Get the data for the second group based on sampled subject matches
		ind2 = c();
		for (s in s2) {
			ind = which(dfIn$Subject == s);
			ind2 = c(ind2,ind);
		}
	
		
		df1 = dfIn[ind1,];
		if (d != 0) {
			df2 = dfCopy[ind2,];
		} else {
			df2 = dfIn[ind2,];
		}
		 
		# Make group labels
		group1 = rep("Group1",length(ind1));
		group2 = rep("Group2",length(ind2));
		
		df1$Group = group1;
		df2$Group = group2;
		
		testDF = rbind(df1,df2);

		# Set plot parameters
		span = 1;	# No smooth since waveforms are smooth already
		tRa = c(-6,6);	# Ranges of t, dof, and p plots
		dfRa = c(5,7);
		pRa = c(-10,0);

		showX = c(T,T,T,T);	# which plots to show scales and legends on
		showY = c(T,T,T,F);
		showL = c(T,F,F,T);
		lPosMean = c(0.7,0.7);	# legend positions
		lPosP = c(0.5,0.25);			
		showPts = T;		# show points on mean plot or no
		lHoriz = T;	# horizontal p legend
				
		gParams = list(0,360,30,"Time Point","Amplitude",span,"TestSineWave",
				tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz)
		
		# Start LAAST;
		p1 = laast(testDF,'X','Y',c('Group1','Group2'),1,'p',gParams);
	
		print(paste('sig = ',exp(min(sig))));
		pvalArr[i] = exp(min(sig));	# presumes that log(pval) was stored
		
		if (pvalArr[i] < alphaN) {
			fp[i] = TRUE;
		} else {
			fp[i] = FALSE;
		}
		ind = which(fp == TRUE);
		rate = length(ind)*100/i;
		print(paste("Current power = ",rate,"%"));

		rateArr[rcnt] <<- rate;
		rcnt <<- rcnt+1;
		
	}	
		
	ind = which(fp == TRUE);
	rate = length(ind)*100/Nrep;
	rateArr[rcnt] <<- rate;
	rcnt <<- rcnt+1;
	if (d == 0) {
		print(paste("False Positive Rate=",rate,"%"));
	} else {
		print(paste("Power=",rate,"%"));
	}

	outString = paste('N1:',N1,',N2:',N2,',d:',d,',Nrep:',Nrep,'Rate:',rate);
	print(outString);
	
	return(outString);

}

# Last edit: May 28, 2019
compareGaussianSamples <- function(Nsamp,Nruns) {
	# This function looks for the false discovery rate of 
	# the Welch T-test on a normal distribution
	# This is a zero-dimensional test which is not used in the paper
	
	d1 = rnorm(1000,0,1);
	
	pval = c();
	for (i in c(1:Nruns) ) {
		s1 = sample(d1,Nsamp,replace=TRUE);
		s2 = sample(d1,Nsamp,replace=TRUE);
		pval[i] = t.test(s1,s2)$p.value;
	}
	
	ind = which(pval < 0.05);
	FDR = length(ind)*100/Nruns;
	print(paste(FDR,"%"));
	
	pvF = data.frame(pval);
	
	p1 = ggplot(pvF,aes(x=pvF))+geom_histogram(alpha=0.3);
	print(p1);

}

# Last edit: May 28, 2019
compareGaussianSamples2 <- function(N1, N2, Npts, d, Nruns, rho) {
	# N1 = number of samples for first sample population
	# N2 = number of samples for second sample population
	# Npts = number of comparisons that would be made along a curve
	# d = effect size or offset between samples
	# Nruns = how many draws to take
	# rho = correlation coefficient 
	# This is a zero-dimensional test which is not used in the paper

	
	# Build normal distribution of 1000 samples with mean of zero and sd of 1.
	d1 = rnorm(10000,0,1);
	
	pval = c();
	for (i in c(1:Nruns) ) {
		s1 = sample(d1,N1,replace=TRUE);
		s2 = sample(d1,N2,replace=TRUE) + d;
		# calcPT <- function(m1,m2,s1,s2,n1a,n2a)
		sd1 = sd(s1);
		sd2 = sd(s2);
		m1 = mean(s1);
		m2 = mean(s2);
		pval[i] = calcPT(m1,m2, sd1,sd2, N1,N2);
		#pval[i] = t.test(s1,s2)$p.value;
	}
	
	alphaAdj = 0.05/( Npts*(1-rho^2)+ rho^2 );
	print(paste("alphaAdj = ",alphaAdj));
	
	ind = which(pval < alphaAdj);
	FDR = length(ind)*100/Nruns;
	
	if (d == 0) {
		print("False discovery rate = ");
	} else { 
		print("Power = ");
	}
	print(paste(FDR,"%"));
	
	pvF = data.frame(pval);
		
}

# Last edit: May 28, 2019
compareNonGaussianSamples <- function(N1, N2, Npts, d, alpha, Nruns, rho) {
	# N1 = number of samples for first sample population
	# N2 = number of samples for second sample population
	# Npts = number of comparisons that would be made along a curve
	# d = effect size or offset between samples
	# alpha = a skewness parameter.  If alpha=0, the distribution should be normal
	#		As alpha increases, skewness does also.
	# Nruns = how many draws to take
	# rho = correlation coefficient 
	# This is a zero-dimensional test which is not used in the paper

	
	# Build normal distribution of 10000 samples with mean of zero and sd of 1.
	#d1 = rnorm(10000,0,1);
	d1 = randn_skew_fast(10000, alpha, loc=0.0, scale=1.0);
		
	cnt = 0;	# n times diff in means equals or exceeds d.
	
	pval = c();
	
	for (i in c(1:Nruns) ) {
		s1 = sample(d1,N1,replace=TRUE);
		s2 = sample(d1,N2,replace=TRUE) + d;
		# calcPT <- function(m1,m2,s1,s2,n1a,n2a)
		sd1 = sd(s1);
		sd2 = sd(s2);
		m1 = mean(s1);
		m2 = mean(s2);
		
		t = calcT(m1,m2, sd1,sd2, N1,N2);	

		if (d != 0) {
			if (t > tcrit) cnt = cnt+1;
		}
		
		pval[i] = calcPT(m1,m2, sd1,sd2, N1,N2);
		#pval[i] = t.test(s1,s2)$p.value;
	}
	
	alphaAdj = 0.05/( Npts*(1-rho^2)+ rho^2 );
	print(paste("alphaAdj = ",alphaAdj));
	
	ind = which(pval < alphaAdj);
	FDR = length(ind)*100/Nruns;
	
	if (d == 0) {
		print("False discovery rate = ");
		print(paste(FDR,"%"));
	} else { 
		print("Power (t-test)= ");
		print(paste(FDR,"%"));
	}
	
	pvF = data.frame(pval);
		
}



# Last edit: May 23, 2019
randn_skew_fast <- function(N, alpha=0.0, loc=0.0, scale=1.0) {
	# Generate a skew normal distribution
	#	N = the number of samples
	#	alpha (this has nothing to do with significance, but with skew - in the text it is given as gamma) = extent of skew
	#	loc = shift of distribution from zero
	# 	scale = scale factor 
	#
	# Based on: 
	# https://stackoverflow.com/questions/36200913/generate-n-random-numbers-from-a-skew-normal-distribution-using-numpy
	#
    sigma = alpha / sqrt(1.0 + alpha^2) 
    u0 = rnorm(N,0,1);	# mean of zero, sd of 1
    v = rnorm(N,0,1);
    
    u1 = (sigma*u0 + sqrt(1.0 - sigma^2)*v) * scale
    ind = which(u0 < 0);
    u1[ind] = -u1[ind]
    #u1[u0 < 0] *= -1
    u1 = u1 + loc
    return (u1);
}

plotRandSkew <- function(N, alpha, loc, scale) {
	u1 = randn_skew_fast(N, alpha, loc, scale);
	
	hist(u1,breaks=50,col='blue');
}
