#	LAASTv2.R
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

# See buildFiguresFinal.R for sample usage.
# 

suppressPlotting = FALSE
plot2Screen <<- TRUE;	# Global, if false, will plot to a file
alphaN <<- c();
alphaH <<- c();
alphaHB <<- c();
	
laast <- function(allData,varX,varY,groups,binSize,type,gParams) {
	# Outputs
	#
	#	Returns a list of 5 ggplot plot objects including: 
	#		mean, t- or f-values, dof, p, and power
	#
	# Inputs
	# 
	#	allData = data frame with time in varX column, and variable varY of interest.
	#			NOTE: a column named 'Group' must exist
	#	varX = generally time variable of allData
	# 	varY = variable of interest of allData
	#	groups = vector containing strings labeling 2 groups found in allData.
	#			NOTE: at the current time LAAST can only deal with two groups at a time.
	#	binSize = related to how frequently data should be sampled.  In a dataset ranging (in x)
	#			from 1-100%, a value of 1 means that sampling will occur at each percentage
	#			of time.  If the binSize is too small to contain a minimum number of points, 
	#			LAAST will not work.
	#	type = type of test conducted and displayed:
	#			'p' = (Default) Welch's t-test showing both t and log p-values
	#			'pp' = Pooled-variance t-test showing both t and log p-values
	#			'f' = F-test showing both F and log p-values
	#			'tp' ????
	#	gParams = list of graphical parameters for plotting
	#		gParams = list(xmin,xmax,deltaX,xtit,ytit,span,fL,
	#			tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts);
	#
	#		xmin = minimum x to display on graph
	#		xmax = maximum x to display on graph
	#		deltaX = increment to display on x-axis of graph
	#		xtit = x-axis title
	#		ytit = y-axis title
	#		span = span for smoothing both in displaying the graph AND in conducting LAAST
	#		fL = file leaf for saving the plot (eg: everything up to the . in 'fileleaf.jpg')
	#		tRa = vector containing the min and max y values to display on t-value or f-value plot
	#		dfRa = vector containing the min and max y values to display on the DOF plot
	#		pRa = vector containing the min and max y values to display on the p-value plot
	#			NOTE that tRa, dfRa, and pRa are needed when multiple panels should plot the
	#			same graph ranges.  If undefined, plotting will try to use default ranges.
	#		showX = vector containing 4 boolean values which indicate whether or not to 
	#			show titles and numbers on the x-axis of the meanPlot, t-plot, dof-plot, and p-plot.
	#		showY = similar vector, but for the y-axis.
	#		showL = similar vector, but for the legend.
	#			NOTE: the show? variables are for multipanel plotting so that scales 
	#			do not end up being redundant
	#		lPosMean = vector giving x and y coords of where the legend for the meanPlot should show up.
	#		lPosP = vector as with lPosMean, but for the p-value plot
	#		showPts = boolean indicating whether or not to show points on the meanPlot
	#
		
	#suppressPlotting is a global

	if (suppressPlotting == F) {	
	
		library(ggplot2);
		library(gridExtra);	# for multipanel plotting with ggplot2
		library(egg);	# Also for multipanel plotting with ggplot2 

	}

	p1 = '';
	p2 = '';
	p3 = '';
	p4 = '';
	p5 = '';
	p6 = '';
	
	# Define globals (aids in limiting number of parameters passed to lower level functions
	
	a1Cnt <<- 1;	# count of alpha-level values from LAAST method
	a2Cnt <<- 1;	# count of alpha-level values from Holm-Bonferroni method
	
	alphaN <<- c();
	alphaH <<- c();
	alphaHB <<- c();
	
	# Plot the mean
	
	if (suppressPlotting == F) {	
		p1 = meanPlot(allData, varX, varY, groups[1],groups[2],gParams,binSize);
	}
	#print('test 1');
	
	if (type == 'p' | type == 't') {
		type1 = 't';
	} else if (type == 'f') {
		type1 = 'f';
	} else {
		type1 = 'tp'
	}

	# Extract some graphical parameters
	xmin = unlist(gParams[1]);
	xmax = unlist(gParams[2]);
	deltaX = unlist(gParams[3]);
	xtit = unlist(gParams[4]);
	ytit = unlist(gParams[5]);
	br = seq(xmin,xmax,deltaX);
	
	tRa = unlist(gParams[8]);
	dfRa = unlist(gParams[9]);
	pRa = unlist(gParams[10]);
	showX = unlist(gParams[11]);
	showY = unlist(gParams[12]);
	showL = unlist(gParams[13]);
	lPosP = unlist(gParams[15]);
	showPts = unlist(gParams[16]);
	lHoriz = unlist(gParams[17]);
	
	# Calc t first and then later p?  Inefficient
	if (suppressPlotting == F) {	
	
		paramsA = smoothAndFit(allData,varX,varY,groups[1],groups[2],binSize,type1,gParams);	# smoothed t-test

		dfAll = getParams(paramsA,type);

		# where dfX = data.frame(times,test,samples,sig,degF);
		dfAll$samples = as.factor(dfAll$samples);
		
		# Plot the degrees of freedom
		p2 = dofPlot(dfAll,gParams);
		p6 = powerPlot(dfAll,gParams);
	}

	# Check ranges.  If tRa != c(0,0), use it.
	# What's tp?  find it.
	if (type == 'p' | type == 'pp' | type == 'tp') {	
		int = log(0.05);
		yra = c(min(sig),max(sig));
		if (max(tRa)-min(tRa) != 0) {
			yra = tRa;
		}
		ytit = 't-values';
	} else {
		int = 0;
		yra = c(0,max(sig));
		ytit = 'F-values';
	}

	if (suppressPlotting == F) {	

		# Set up theme
		th = theme();
		if (showY[3]) {
			if (showX[3]) {
				theme3 = th;
			} else {
				theme3 = th + theme(axis.title.x=element_blank(), 
					axis.text.x=element_blank());		
			}
		} else {
			if (showX[3]) {
				theme3 = th + theme(axis.title.y=element_blank(), 
					axis.text.y=element_blank()) 					
			} else {
				theme3 = th + theme(axis.title.x=element_blank(), 
					axis.text.x=element_blank(), 
					axis.title.y=element_blank(), 
					axis.text.y=element_blank()) 		
			}
		}
			
		# Plot the t-values
		p3 = ggplot(dfAll)+geom_line(aes(x=times,y=sig), color="#10BE44",size=0.75,linetype="longdash" )+
			#geom_abline(intercept=int,slope=0)+
			ylab(ytit)+
			scale_x_continuous(breaks=br,name=xtit)+
			ylim(yra[1],yra[2])+theme3;
	}		
	# Now get p-values
	if (type == 'p' | type == 't' | type == 'f') {
		type1 = 'p'; # Use Welch
	} else {
		type1 = 'pp'
	}

	paramsA = smoothAndFit(allData,varX,varY,groups[1],groups[2],binSize,type,gParams);	# smoothed t-test

	dfAll = getParams(paramsA,type);
	sig <<- dfAll$sig
	dfAll$samples = as.factor(dfAll$samples);		
	
	CIlower = dfAll$CIlower;	#dfA vs dfAll?
	CIupper = dfAll$CIupper;
	times = dfAll$times;
	deltaM = dfAll$deltaM;
	
	CI = data.frame(times,deltaM,CIlower,CIupper);

	int = log(0.05);	# log(0.05);
	print(paste("**** LAAST alphas",a1Cnt));
	print(alphaN);
	print(paste("**** Holm alphas",a2Cnt));
	print(alphaH);
	print(paste("**** Hochberg alphas",a2Cnt));
	print(alphaHB);

	# End points of horizontal lines showing adjusted alpha levels
	int1 = log(alphaN[1]);	# LAAST alpha values
	int2 = log(alphaH[1]);	# Holm-Bonferroni alpha values
	int3 = log(alphaHB[1]);	# Hochberg alpha values

	# Set up labeling for alpha value lines
	# 	Note that the order here doesn't much matter since R will automatically
	#	put legends in alphabetical order unless specifically told not to.
	#   So one just needs to get the order consistent between vectors that
	# 	comprise the dataframe legDF.
	#
	# Default values if no RFT shown
	ints = c(int1,int2,int3);
	int1vec = rep(int1,2);	# only need end points
	int2vec = rep(int2,2);
	int3vec = rep(int3,2);
	
	lab1vec = rep("LAAST",2);
	lab2vec = rep("Holm.",2);
	lab3vec = rep("Hoch.",2);
	
	lts = c("dotted", "longdash","dashed");	
	cols = c("purple","red","#20C250");	
	
	intVec = c(int1vec,int2vec,int3vec);
	labVec = c(lab1vec,lab2vec,lab3vec);
	
	xMin = min(dfAll$times);
	xMax = max(dfAll$times);	
	xV = c(xMin,xMax,xMin,xMax,xMin,xMax);
	
	
	legDF = data.frame(xV,intVec,labVec);
	names(legDF)[3] = "Alpha";
	
	# New ytit for next graph
	ytit = 'log(p-value)';
		
	yra = c(min(sig),0);
	if (max(pRa)-min(pRa) != 0) {
		yra = pRa;
	}
	
	#legend.position="none"
	if (suppressPlotting == F) {		
		# Set up theme
		if (lHoriz == T) {
			th = theme(legend.direction="horizontal",legend.position=lPosP,legend.background=element_blank());
		} else {
			th = theme(legend.direction="vertical",legend.position=lPosP,legend.background=element_blank());		
		}
		if (showY[4]) {
			if (showL[4]) {
				theme4 = th;
			} else {
				theme4 = th + theme(legend.position="none");
			}
		} else {
			if (showL[4]) {
				theme4 = th + theme(axis.title.y=element_blank(), 
					axis.text.y=element_blank());	
			} else {
				theme4 = th + theme(legend.position="none",
					axis.title.y=element_blank(), 
					axis.text.y=element_blank());	
			}
		}
		
	
		# p-value plots
		p4 = ggplot(dfAll)+
				geom_line(aes(x=times,y=sig), size=1.25  )+
				geom_line(data=legDF,aes(x=xV,y=intVec,color=Alpha,linetype=Alpha))+
				scale_x_continuous(breaks=br,name=xtit)+ylab(ytit)+
				scale_linetype_manual(values=lts)+
				scale_color_manual(values=cols)+
				ylim(yra[1],yra[2])+theme4;
	
		# Confidence Interval Plots
		p5 = ggplot(CI) + geom_line(aes(x=times,y=deltaM),size=1.25,color="black") +
			geom_line(aes(x=times,y=CIlower),size=1.25,color="red") +
			geom_line(aes(x=times,y=CIupper), size=1.25, color="red") + 
			geom_abline(intercept=0,slope=0,linetype=2,colour="black")+
			ylab('95% CI')+
			scale_x_continuous(breaks=br,name=xtit);
		
		g1 <- ggarrange(p1, p3, p2, p4, nrow = 2);	# from egg - better for 
			# Controlling plot widths	

		if (plot2Screen) {
			print(g1);
		} else {
			
			fn = unlist(gParams[7]);
			fN = paste(fn,".tiff",sep="");
			print(paste("Plotting ",fN,sep=""));
			#tiff(filename = fN, pointsize =12, bg = "white", res = 500)
			ggsave(filename=fN,plot=g1,width=8,height=6,dpi=500);  # units="in",
			#dev.off();
		}
	}
	plotsOut <- list(p1,p2,p3,p4,p5,p6);
	
	return (plotsOut);
}


# Organize parameters from smoothAndFit for easier plotting with ggplot2
getParams <- function(params,type) {
	times = params$times;
	sig = params$pvals; 
	degF = params$degF;
	N = length(times);
	samples = c();
	samples[1:N] = N;
	test = c();
	test[1:N] = type;
	m1 = params$m1;
	m2 = params$m2;
	deltaM = m1-m2;
	CIlower = params$CIlower;
	CIupper = params$CIupper;
	pwr = params$pwr;
	
	dF = data.frame(times,test,samples,sig,degF,deltaM,CIlower,CIupper,pwr);
	return (dF);
		
}

# Plot effective number in groups vs time
nPlot <- function(n1,n2,times,groups) {
	Times = c(times,times);
	N = c(n1,n2);
	l1 = length(n1);
	l2 = length(n2);

	group1 = c();
	group1[1:l1] = groups[1];
	group2 = c();
	group2[1:l2] = groups[2];
	
	Group = c(group1,group2);

	nData = data.frame(Times,N,Group);
	
	p2 = ggplot(nData,aes(x=Times,y=N,color=Group)) + geom_line() + xlab("% Cycle")+
		ylab("Numbers in Group");
		
	return (p2);
	
}

# Plot degrees of freedom vs time
dofPlot <- function(dfAll,gParams) {
	xtit = unlist(gParams[4]);

	xmin = unlist(gParams[1]);
	xmax = unlist(gParams[2]);
	deltaX = unlist(gParams[3]);
	xtit = unlist(gParams[4]);
	dfRa = unlist(gParams[9]);
	showX = unlist(gParams[11]);
	showY = unlist(gParams[12]);
	
	#ytit = unlist(gParams[5]);
	#print(gParams);
	br = seq(xmin,xmax,deltaX);

	th = theme();
	if (showY[2]) {
		if (showX[2]) {
			theme2 = th;
		} else {
			theme2 = th + theme(axis.title.x=element_blank(), 
				axis.text.x=element_blank()); 
		}
	} else {
		if (showX[2]) {
			theme2 = th+ theme(axis.title.y=element_blank(), 
				axis.text.y=element_blank());
		} else {
			theme2 = th + theme(axis.title.y=element_blank(), 
				axis.text.y=element_blank(),
				axis.title.x=element_blank(), 
				axis.text.x=element_blank()); 
		}		
	}

	if (max(dfRa)-min(dfRa) != 0) {
		dfRa = dfRa;	
	} else {
		dfRa = c(min(dfAll$degF),max(dfAll$degF));
	}
	
	# Check to see if range was manually set
	if (max(dfRa)-min(dfRa) != 0) {	# yes it was
		p2 = ggplot(dfAll,aes(x=times,y=degF)) + geom_line(size=1.25,color="black") + 
			scale_x_continuous(breaks=br,name=xtit)+ylab("DOF")+
			ylim(dfRa[1],dfRa[2])+ theme2;
		
	} else {	# No it wasn't
		p2 = ggplot(dfAll,aes(x=times,y=degF)) + geom_line(size=1.25,color="red") + 
			scale_x_continuous(breaks=br,name=xtit)+ylab("DOF")+theme2;
	}
	
	return (p2);
}

# Plot power vs time - not used in paper
powerPlot <- function(dfAll,gParams) {
	xtit = unlist(gParams[4]);

	xmin = unlist(gParams[1]);
	xmax = unlist(gParams[2]);
	deltaX = unlist(gParams[3]);
	xtit = unlist(gParams[4]);
	#ytit = unlist(gParams[5]);
	#print(gParams);
	br = seq(xmin,xmax,deltaX);
	
	p6 = ggplot(dfAll,aes(x=times,y=pwr,color=samples)) + geom_line(size=1.25) + 
		scale_x_continuous(breaks=br,name=xtit)+ylab("Power");
		#xlab(xtit)+ylab("Power");
	
	return (p6);
}

# Plot mean and SD of groups vs time
meanPlot <- function(dataS, varX, varY, group1, group2, gParams,binSize0) {
	
		
	# Index which group is being modeled
	i1 = which(dataS$Group == group1 | dataS$Group == group2); # Model specific condition (Group)	

	data1 = dataS[i1,];	# subset data of interest from dataS

	# Extract x and y variables from subsetted data
	x= data1[varX][,1];
	y= data1[varY][,1];
	
	x1 = x;	# Points for plotting
	y1 = y;
	
	group = data1["Group"][,1];
	
	# Put it into a data frame
	df0 = data.frame(x,y,group);
	
	xmin = unlist(gParams[1]);
	xmax = unlist(gParams[2]);
	deltaX = unlist(gParams[3]);
	xtit = unlist(gParams[4]);	
	ytit = unlist(gParams[5]);

	showX = unlist(gParams[11]);
	showY = unlist(gParams[12]);
	legPosMean = unlist(gParams[14]);
	showPts = unlist(gParams[16]);
	
	uG = unique(data1$Group);
	
	##

	mod1 <- modelSeries(dataS, varX, varY, group1, binSize0,gParams)	
	#print(mod1);
	m1 = mod1$m1;
	se1a = mod1$sd1;
	n1a = mod1$n1a;
	#g1 = rep('Group1',length(m1));
	g1 = rep(group1,length(m1));
	
	mod2 <- modelSeries(dataS, varX, varY, group2, binSize0,gParams)	
	m2 = mod2$m1;
	se2a = mod2$sd1;
	n2a = mod2$n1;
	
	#g2 = rep('Group2',length(m2));
	g2 = rep(group2,length(m2));
	
	br = seq(xmin,xmax,deltaX);
	nsamp = (xmax-xmin)/binSize0;
	xArr = c(1:nsamp)*binSize0+ xmin;	# Well behaved x values for model

	y = c(m1,m2);
	ymin = c(m1-se1a,m2-se2a);
	ymax = c(m1+se1a,m2+se2a);
	sdl = c(se1a,se2a);
	group = c(g1,g2);
	x = c(xArr,xArr);
	
	df1 = data.frame(x,y,ymin,ymax,sdl,group);
	
	th = theme(legend.position=legPosMean,legend.background=element_blank());

	if (showX[1]) {
		if (showY[1]) {
			th = th;
		} else {
			th =  th + theme(axis.title.y=element_blank(), 
				axis.text.y=element_blank());		
		}
	} else {
		if (showY[1]) {
			th = th + theme(axis.title.x=element_blank(), 
				axis.text.x=element_blank());
		} else {
			th =  th + theme(axis.title.y=element_blank(), 
				axis.text.y=element_blank(),
				axis.title.x=element_blank(),
				axis.text.x=element_blank());		
		}		
	}


	if (showPts) {
		p1 = ggplot(df1,aes(x=x,y=y,color=group,fill=group,linetype=group))+
			#geom_smooth(size=1)+ylab(ytit)+  # this line does se rather than sd
			#stat_summary(fun.data="mean_sd",geom="smooth",size=1)+ylab(ytit)+
			#stat_smooth(geom="smooth",se=T)+
			geom_ribbon(aes(ymin=ymin, ymax=ymax, x=x), alpha = 0.3)+
			ylab(ytit)+
			geom_line()+
			geom_point(data=df0,aes(x=x,y=y,color=group),size=0.4,alpha=0.4)+
			scale_x_continuous(breaks=br,name=xtit)+th;
			#theme(legend.position=legPosMean);

	} else {
		p1 = ggplot(df1,aes(x=x,y=y,color=group,fill=group,linetype=group))+
			#geom_smooth(size=1)+ylab(ytit)+  # this line does se rather than sd
			#stat_summary(fun.data="mean_sd",geom="smooth",size=1)+ylab(ytit)+
			#stat_smooth(geom="smooth",se=T)+
			geom_ribbon(aes(ymin=ymin, ymax=ymax, x=x), alpha = 0.3)+
			ylab(ytit)+
			geom_line()+
			scale_x_continuous(breaks=br,name=xtit)+th;
			#theme(legend.position=c(0.1,0.9));
	}
	return (p1);
}

smoothAndFit <- function(dataS, varX,varY,group1,group2,binSize0,type,gParams) {
	
	# dataS is a data frame containing columnes of data with names varA, varB...
	# varX is a string of the x axis variable name from dataS
	# varY is a string of the y axis variable name from dataS
	# group1,group2 are the grouping variables for either the "Group" column or the "Gender" 
	# 	column of dataS
	
	# binSize0 is the bin size along the x axis.  This determines how much data
	# 	are grouped together for analysis.
	
	# Load needed R libraries for this analysis.  

	
	library(fields);

	x= dataS[varX][,1];
	#xmax = 100;
	#xmin = 1;
	xmin = unlist(gParams[1]);
	xmax = unlist(gParams[2]);
	nsamp = (xmax-xmin)/binSize0;
	times = c(1:nsamp)*binSize0+ xmin;	
	
	
	mod1 <- modelSeries(dataS, varX, varY, group1, binSize0,gParams)	
	m1 = mod1$m1;
	se1a = mod1$sd1;
	n1a = mod1$n1;
	#print('************* n1a **************')
	#print(n1a);
	
	mod2 <- modelSeries(dataS, varX, varY, group2, binSize0,gParams)	
	m2 = mod2$m1;
	se2a = mod2$sd1;
	n2a = mod2$n1;
	#print('************* n2a **************')
	#print(n2a);
	
	#getPower2n <- function(m1, m2, s1, s2, n1, n2) 
	#print("Getting power");
	pwr = getPower2n(m1, m2, se1a, se2a, n1a, n2a, type);
	
	# Use calcPT for p values and
	# calcT for t values	
	alpha <- 0.05;	# Global - changed to local May 31, 2018
	
	if (type == 'p') {
		pvals = calcPT(m1,m2,se1a,se2a,n1a,n2a);
		alphaH[a2Cnt] <<- calcPCritHolm(pvals);
		alphaHB[a2Cnt] <<- calcPCritHochberg(pvals);
		a2Cnt <<- a2Cnt+1;
		
		print(paste('Pcrit Holm = ',alphaH[a2Cnt-1]));

		pvals = log(pvals);	# take log for contrast
		
		alphaAdj = TRUE;
		if (alphaAdj) {
			t = calcT(m1,m2,se1a,se2a,n1a,n2a);	
			t = t[!is.na(t)];
			
			tL = t[1:length(t)-1];
			tR = t[2:length(t)];
			rho = cor(tL,tR);
			
			# Use something like Bonferroni correction, but adjust
			# according correlation coefficient.  If rho == 0, 
			# use full Bonferroni.  If rho == 1, alpha == 0.05
			alphaN[a1Cnt] <<- 0.05/( length(t)*(1-rho^2)+ rho^2 );
			a1Cnt <<- a1Cnt+1;
			
			print(paste('Alpha = ',alpha,", rho = ",rho,", length(t) = ",length(t)) );
		}
		
		
	} else if (type == 't') {
		pvals = calcT(m1,m2,se1a,se2a,n1a,n2a);
	} else if (type == 'pp') {
		# p-pooled sd
		#n1a = unlist(gParams[8]);
		#n2a = unlist(gParams[9]);
		
		sp = getPooledSD(se1a,se2a,n1a,n2a);
		pvals = calcPTpooled(m1,m2,sp,n1a,n2a);

		alphaHB[a2Cnt] <<- calcPCritHochberg(pvals);
		alphaH[a2Cnt] <<- calcPCritHolm(pvals);
		a2Cnt <<- a2Cnt+1;
		
		#print(paste('Pcrit Holm = ',pc));

		pvals = log(pvals);	# take log for contrast
		
		alphaAdj = TRUE;
		if (alphaAdj) {
			t = calcPTpooled(m1,m2,sp,n1a,n2a);	
			t = t[!is.na(t)];
			
			tL = t[1:length(t)-1];
			tR = t[2:length(t)];
			rho = cor(tL,tR);
			
			# Use something like Bonferroni correction, but adjust
			# according correlation coefficient.  If rho == 0, 
			# use full Bonferroni.  If rho == 1, alpha == 0.05
			alphaN[a1Cnt] <<- 0.05/( length(t)*(1-rho^2)+ rho^2 );
			a1Cnt <<- a1Cnt+1;
			
			print(paste('Alpha = ',alpha,", rho = ",rho,", length(t) = ",length(t)) );
		}

	} else if (type == 'tp') {
		
		sp = getPooledSD(se1a,se2a,n1a,n2a);
		pvals = calcTpooled(m1,m2,sp,n1a,n2a);
		
	} else {
		pvals = calcT(m1,m2,se1a,se2a,n1a,n2a);
		pvals = pvals*pvals;	# F test is t test squared		
	}

	# recalc DOF to return
	if (type == 'p' | type == 't') {
		degF = calcDF(se1a,se2a,n1a,n2a); # Welch Satterthwaite
	} else {
		degF = round(n1a + n2a - 2);			# Pooled, no SD dependence
	}
	
	CI = calcWCI(0.05,m1,m2,se1a,se2a,n1a,n2a)
	CIlower = CI$CIlower;
	CIupper = CI$CIupper;
	
	pDF = data.frame(times,pvals,m1,m2,se1a,se2a,n1a,n2a,degF,CIlower,CIupper,pwr);

	return(pDF);
}


modelSeriesXY <- function(dfXY, binSize0,gParams) {
	# Differs from modelSeries in that the x, y variables have been selected
	# 	as well as the group
	
	# dataS is a data frame containing columnes of data with names varA, varB...
	# varX is a string of the x axis variable name from dataS
	# varY is a string of the y axis variable name from dataS
	# group1 is the grouping variable for either the "Group" column or the "Gender" 
	# 	column of dataS
	
	# binSize0 is the bin size along the x axis.  This determines how much data
	# 	are grouped together for analysis.
	
	# Load needed R libraries for this analysis.  
	library(fields);
	library(ggplot2);
	library(msir);

	# span0 is a parameter used in loess for smoothing.  Larger spans have smoother results
	# Default is 0.75
	#span0 = 0.75
	#span0 = 0.10; # for Pataky
	span0 = unlist(gParams[6]);
	
	# We are apriori assuming that we're dealing with gait data that goes from 1-100%
	# If varX is well behaved, one could pull these values from it instead	
	#xmax = 100;
	#xmin = 1;
	xmax = unlist(gParams[2]);
	xmin = unlist(gParams[1]);
	nsamp = (xmax-xmin)/binSize0;
	
	# Number of bins
	# 	used to start at 0
	xArr = c(1:nsamp)*binSize0+ xmin;	# Well behaved x values for model

	# Smooth/model it using loess
	
	#print('test1');
	sp0 = loess(y~x,dfXY,span=span0);
	#print('test2');
	# Predict central values based on xArr
	sp0P = predict(sp0,xArr,se=TRUE);

	m1 = sp0P$fit;	# mean of model to data at values xArr
	
	# Get standard deviation of the data
	x = dfXY$x;
	y = dfXY$y;
	
	lsd = loess.sd(x,y,nsigma=1,span=span0);

	y1 = lsd$sd;
	x1 = lsd$x;
	#plot(x1,y1);
		
	# Put sd values into a data frame	
	df1 = data.frame(x1,y1);
	# Smooth/model sd values using loess
	sp1 = loess(y1~x1,df1,span=span0);
	# Predict sd values based on xArr
	sp1P = predict(sp1,xArr,se=TRUE);
	sd1 = sp1P$fit;  # modeled standard deviation of data at values xArr

	# Get CI envelope and estimate
	
	# Estimate number of data points at each sample based on density of x
	nsamp = (xmax-xmin)/binSize0
	#nsamp = (100-0)/binSize0;	# Number of samples
	n1 = nPerSamp(x,nsamp,xmin,xmax);
	
	#n1 = nPerBin(x,1,100,binSize0);

	ymin = min(y);
	ymax = max(y);
	
	doPlot = FALSE;
	if (doPlot) {
		plot(x,y,col="blue",ylim=c(ymin,ymax),main=group1);
		lines(x,lsd$y,col="blue");
		lines(x,lsd$upper,col="red");
		lines(x,lsd$lower,col="red");
	}
		
	#print(length(n1));
	#print(length(m1));
	#print(length(sd1));
	modelOut <- data.frame(m1,n1,sd1);
	return (modelOut);
	
}

modelSeries <- function(dataS, varX, varY, group1, binSize0,gParams) {
	
	# dataS is a data frame containing columns of data with names varA, varB...
	# varX is a string of the x axis variable name from dataS
	# varY is a string of the y axis variable name from dataS
	# group1 is the grouping variable for either the "Group" column or the "Gender" 
	# 	column of dataS
	
	# binSize0 is the bin size along the x axis.  This determines how much data
	# 	are grouped together for analysis.
	
	# Load needed R libraries for this analysis.  
	library(fields);
	library(ggplot2);
	library(msir);

	# span0 is a parameter used in loess for smoothing.  Larger spans have smoother results
	# Default is 0.75
	#span0 = 0.75
	#span0 = 0.10; # for Pataky
	span0 = unlist(gParams[6]);
	
	# We are apriori assuming that we're dealing with gait data that goes from 1-100%
	# If varX is well behaved, one could pull these values from it instead
	x= dataS[varX][,1];
	
	#xmax = 100;
	#xmin = 1;
	xmax = unlist(gParams[2]);
	xmin = unlist(gParams[1]);
	nsamp = (xmax-xmin)/binSize0;
	
	# Number of bins
	# 	used to start at 0
	xArr = c(1:nsamp)*binSize0+ xmin;	# Well behaved x values for model

	# Index which group is being modeled based on either Gender or Group
	if (group1 == "Male" | group1 == "Female") {
		i1 = which(dataS$Gender == group1);	# Model specific gender
	} else {
		i1 = which(dataS$Group == group1); # Model specific condition (Group)
	}

	data1 = dataS[i1,];	# subset data of interest from dataS

	# Extract x and y variables from subsetted data
	x= data1[varX][,1];
	y= data1[varY][,1];

	#print('printing df');
	#print(x);
	#print(y);
	# Put it into a data frame
	df0 = data.frame(x,y);
	# Smooth/model it using loess
	
	#print(df0);
	#print('test1');
	sp0 = loess(y~x,df0,span=span0);
	#print('test2');
	# Predict central values based on xArr
	sp0P = predict(sp0,xArr,se=TRUE);

	m1 = sp0P$fit;	# mean of model to data at values xArr
	
	# Get standard deviation of the data

	lsd = loess.sd(x,y,nsigma=1,span=span0);

	y1 = lsd$sd;
	x1 = lsd$x;
	#plot(x1,y1);
		
	# Put sd values into a data frame	
	df1 = data.frame(x1,y1);
	# Smooth/model sd values using loess
	sp1 = loess(y1~x1,df1,span=span0);
	# Predict sd values based on xArr
	sp1P = predict(sp1,xArr,se=TRUE);
	sd1 = sp1P$fit;  # modeled standard deviation of data at values xArr

	# Get CI envelope and estimate
	#print(group1);
	# Estimate number of data points at each sample based on density of x
	nsamp = (xmax-xmin)/binSize0
	#nsamp = (100-0)/binSize0;	# Number of samples
	n1 = nPerSamp(x,nsamp,xmin,xmax);
	
	#n1 = nPerBin(x,1,100,binSize0);

	ymin = min(y);
	ymax = max(y);
	
	doPlot = FALSE;
	if (doPlot) {
		plot(x,y,col="blue",ylim=c(ymin,ymax),main=group1);
		lines(x,lsd$y,col="blue");
		lines(x,lsd$upper,col="red");
		lines(x,lsd$lower,col="red");
	}
		
	#print(length(n1));
	#print(length(m1));
	#print(length(sd1));
	modelOut <- data.frame(m1,n1,sd1);
	return (modelOut);
	
}

getPower2n <- function(m1, m2, s1, s2, n1, n2, type) {
	#pwr.t2n.test(n1 = , n2= , d = , sig.level =, power = )
	library(pwr);
	#print("Diff in means");
	#print(abs(m1-m2));
	
	d = getEffectSize(m1, m2, s1, s2, n1, n2);
	#print("Effect sizes");
	#print(d);
	pow = pwr.t2n.test(n1 = n1, n2 = n2, d = d, sig.level = 0.05);
	
	powOut = pow$power; 
	
	#print("Powers");
	#print(powOut);
	return (powOut);
	
}

getEffectSize <- function(m1,m2,s1,s2,n1,n2) {
	
	sP = getPooledSD(s1,s2,n1,n2);
	
	#print("Pooled SD");
	#print(sP);

	# Cohen's d
	d = abs(m1-m2)/sP;
	
	return (d);
}

getPooledSD <- function(s1,s2,n1,n2) {
	
	s2 = ( (n1-1)*s1*s1 + (n2-1)*s2*s2 ) / (n1 + n2 - 2);
	
	s = sqrt(s2);
	
	return (s)
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


nPerBin <- function(timePts,N,xmin,xmax) {
	# Estimate the number of points in each BIN along the time domain.
	# This function uses a density function to model the frequencies
	# of data across all time points.  This density estimate is then
	# resampled at N intervals and scaled up to the total number of
	# original data points.  
	
	
	# Extend timePts to avoid edge effects
	timePts = sort(timePts);	# has to be in order or this doesn't work	
	
	b = -timePts;
	dba = min(timePts)-max(b);
	rTPleft = -rev(timePts)+dba-1
	
	dbb = max(timePts)-min(b);
	rTPright = rev(-timePts)+dbb+1
	
	x = c(rTPleft,timePts,rTPright);
	
	Nt = length(timePts);
	
	#print('############### from nPerSamp ######################');

	d = density(x,na.rm=TRUE);
	
	d1x = d$x;
	d1y = d$y
	n = length(timePts);
	d1y = n*d1y/3; # scale to total number of points
	
	ind = which(d1x > xmin & d1x < xmax);
	d1x = d1x[ind];
	d1y = d1y[ind];
	 
	appr = approx(d1x,d1y,n=N);  # Resampling or approximation
	
	d2y = appr$y;
	
	ratio = Nt/N;	# Ratio of total number of timePts to sample number
	d2y = ratio*d2y/2;	# Too large by factor of 2 - Jan 22, 2018
	
	ind = which(d2y < 1);
	d2y[ind] = 1;	# Need at least 1 data point per sample
	
	return (d2y);
	
}


sampleEvery <- function(vec, N) {
	l = length(vec);
	vOut = c();
	dx = l/N;
	cnt = 1;
	for (i in 1:N) {
		vOut[cnt] = vec[i*dx];
		cnt = cnt+1;
	}
	return (vOut);
}

calcPCritHolm <- function(pvals) {
	m = length(pvals);
	pS = sort(pvals);
	pcrit = 0.05/m;
	cnt = 1;
	for (p in pS) {
		if (pS[cnt] < 0.05/m) {
			
			pcrit = 0.05/m;
			m = m-1;
		}
		#print(paste('pcrit = ',pcrit));
		cnt = cnt+1;
	}
	return (pcrit);
}

calcPCritHochberg <- function(pvals) {
	m = length(pvals);
	pS = sort(pvals);
	pcrit = 0.05/m;
	k = 1;
	for (p in pS) {
		if (pS[k] < 0.05*k/m) {
			
			pcrit = 0.05*k/m;
			#m = m-1;
		}
		#print(paste('pcrit = ',pcrit));
		k = k+1;
	}
	return (pcrit);	
}

getDFstd <- function(allData,group1,group2) {
	ind1 = which(allData$Group == group1);
	ind2 = which(allData$Group == group2);
	
	n1 = length(ind1);
	n2 = length(ind2);
	
	df1 = calcDFstd(n1,n2);
	
	return(df1);
}

calcDFstd <- function(n1,n2) {
	v1 = n1-1;
	v2 = n2-1;
	df1 = v1+v2;
	return (df1)
}

calcDF <- function(s1,s2,n1,n2) {
	# Using the Welch-Satterwaite equation,  calculate the degrees of freedom
	# necessary for a Welch Two Sample t-test
	
	v1 = n1-1;
	v2 = n2-1;
	
	# Welch-Satterthwaite Equation for degrees of freedom
	df1 = ( (s1^2)/n1 + (s2^2)/n2)^2/( (s1^4)/(n1*n1*v1) + (s2^4)/(n2*n2*v2))

	return (df1);
}

calcPT <- function(m1,m2,s1,s2,n1a,n2a) {
	
	# Calculate p(t) - probability of t-value given means, sds, and numbers in 
	# each distribution of interest.  Uses Welch 2-sample t-test.

	# Actually calculate T
	t = abs(calcT(m1,m2,s1,s2,n1a,n2a));  # changed on Dec 31, 2017	
	
	# Calculate degrees of freedom
	df1 = calcDF(s1,s2,n1a,n2a);
			
	pOut = 2*pt(t,df1,lower=FALSE);
	
	return (pOut);
}

# return t values, not p values
calcT <- function(m1,m2,s1,s2,n1a,n2a) {
	#print('Calling calcT');
	# Calculate t values given means, sds, and numbers in 
	# each distribution of interest.  Uses Welch 2-sample t-test.
	
	# Welch t-test
	t = (m1 - m2) / sqrt( (s1^2)/n1a + (s2^2)/n2a )
	
	sig = sqrt( (s1^2)/n1a + (s2^2)/n2a );
	#### Fix Me
	#t = abs(t);	# two tailed

	return (t);
}

calcPTpooled <- function(m1,m2,s,n1a,n2a) {
	
	# Calculate p(t) - probability of t-value given means, sds, and numbers in 
	# each distribution of interest.  Uses Pooled 2-sample t-test.

	# Actually calculate T
	t = abs(calcTpooled(m1,m2,s,n1a,n2a));  # Changed on Feb 1, 2018	

	# Calculate degrees of freedom
	df1 = n1a + n2a - 2;
	
	#pOut = 2*pt(t,df1,lower=FALSE);
	pOut = 2*pt(t,df1,lower=FALSE);
	
	return (pOut);
}

calcTpooled <- function(m1,m2,s,n1a,n2a) {
	#print('Calling calcTpooled');
	# Calculate t-values given means, pooled sd, and numbers...
	t = (m1 - m2) / (s* sqrt( 1/n1a + 1/n2a ) )
	#t = abs(t);	# two tailed

	return (t);
	
}

calcWCI <- function(alpha,m1,m2,s1,s2,n1a,n2a) {
	# Calculate Welch confidence interval
	
	DOF = calcDF(s1,s2,n1a,n2a);
	sd12 = sqrt( (s1^2)/n1a + (s2^2)/n2a );
	
	p = 1-alpha/2;
	t = qt(p,DOF);
	
	delta = m1-m2

	CIupper = delta+t*sd12;
	CIlower = delta-t*sd12;
	
	#print("begin CI");
	#print(CIupper);
	#print("end CI");
	
	CI = data.frame(CIlower,CIupper);
	
	return (CI);
}

cherryPick <- function(dataS, varY, groups) {
	# groups is an array c();
	subjects = unique(dataS$Subject);

	y= dataS[varY][,1];
	#print(y);
	maxVals1 = c();
	maxVals2 = c();
	cnt1 = 1;
	cnt2 = 1;
	for (s in subjects) {
		ind1 = which(dataS$Subject == s & dataS$Group == groups[1]);
		if (length(ind1) != 0) {
			yi1 = y[ind1];
			maxVals1[cnt1] = max(yi1);	# Pick the maximum;
			cnt1 = cnt1 + 1;
		}
		ind2 = which(dataS$Subject == s & dataS$Group == groups[2]);
		if (length(ind2) != 0) {
			yi2 = y[ind2];
			maxVals2[cnt2] = max(yi2);	# Pick the maximum;
			cnt2 = cnt2 + 1;
		}
	}
	print(maxVals1);
	print(maxVals2);
	print(groups);
	print(t.test(maxVals1,maxVals2));
		
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plotFinalGraph <- function(tVals4plot) {
	intN = qt(0.0168/2,df=42,lower.tail=FALSE);
	#print(intN);
	
	intP = 2.77;	# From Pataky et al (2015), Fig 7c
	
	ggplot(tVals4plot,aes(x=Time,y=tval,color=Statistic))+
		geom_line(size=1)+xlab("%Gait Cycle")+ylab("t-values")+
		geom_abline(intercept=intP,slope=0,linetype=2,colour="blue",alpha=0.3)+
		geom_abline(intercept=intN,slope=0,linetype=2,colour="red",alpha=0.3);
		
		# 2.77 is tcrit for Pataky
		
}

# Store these as globals to evaluate smoothing
rmsArr = c(); # RMS values found in the smoothing iterations
rmsCnt = 1;	  # index of rmsArr and spanArr
spanArr = c();# span values found in the smoothing iterations

findSpan <- function(x,y) {
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
