# buildFiguresFinal.R

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


##### SUMMARY ######
	##	Note that this code depends on LAASTv2.R which should be distributed alongside this file.
	##

	##### Functions Used in Building Simulated Data or Reading in Real Data #####
	
	# data1a - 		buildData1a() - 			simulated data used in Figure 1
	# patakyData -	readPatakyData() - 			real data used in Figure 2

	##### Functions Used in Building Figures for Paper #####

	# Graphical Abstract - graphicalAbstract()
	# Figure 1 - 	buildFigure1()
	# Figure 2 - 	buildFigure2()

	###### Functions Used in Building Figures & Tables in the Appendix #####
	
	# Figure B.1  -	runBinDependenceBumpy() -	produce graph showing how adjusted alpha level varies
	# Figure B.2 - 	alphaPatterns()
	# Figure B.3  -	runBinDependenceBumpy() - 	produce graph showing example of data and regions of significance.
	# Table C.1  -	See file "buildDistributions.R"
	# Table C.2  -	See file "buildDistributions.R"
	
##### END SUMMARY #####

# This function is used in generating data for Figure 1
buildData1a <- function(f, phaseOffset) {
	# Builds two random sin based sequence GROUPS for comparison
	# Both have 360 points per curve per group so RFT can LAAST can
	#  be used for comparison.
	#
	# f is a scale factor for the 2nd sinusoid
	# phaseOffset = boolean, if true, use gaussian phase angle offset

	Xrad = (pi/180)*c(1:360);
	Y = sin(Xrad);
	
	X1 = c();
	X2 = c();
	Y1 = c();
	Y2 = c();
	
	Yout = c();

	for (i in 1:20) {
		#rad = 10*runif(1)*pi/180;	# Random phase angle scatter
		rad = 20*rnorm(1,0,1)*pi/180;	# Gaussian phase angle scatter
		if (phaseOffset) {
			rad = 0;
		}
		Y = sin(Xrad + rad)+rnorm(1,0,1);	# Gaussian offset
		Y1 = c(Y1,Y);
		Y = f*sin(Xrad + rad)+rnorm(1,0,1);
		Y2 = c(Y2,Y);
		
		X = c(1:360);
		X1 = c(X1,X);
		X2 = c(X2,X);
	}
	l = length(X1);
	
	name1 = rep("Group1",l);
	name2 = rep("Group2",l);
	
	X = c(X1,X2);
	Y = c(Y1,Y2);
	Group = c(name1, name2);
	
	testDF = data.frame(X,Y,Group);
	
	testDFa = data.frame(X1,Y1,name1);
	names(testDFa) = c("X","Y","Group");
	if (phaseOffset) {
		fn = paste("dataSineTestPhase",f,"a.txt",sep="");
	} else {
		fn = paste("dataSineTest",f,"a.txt",sep="");
	}
	write.table(testDFa,fn,sep="\t",quote=FALSE,row.names=FALSE);

	testDFb = data.frame(X2,Y2,name2);
	names(testDFb) = c("X","Y","Group");
	if (phaseOffset) {
		fn = paste("dataSineTestPhase",f,"b.txt",sep="");
	} else {
		fn = paste("dataSineTest",f,"b.txt",sep="");
	}
	
	# Write an output file for use by testRFT.py
	write.table(testDFb,fn,sep="\t",quote=FALSE,row.names=FALSE);

	return (testDF);
}


# Last edit: June 3, 2019
readPatakyT <- function() {
	# Read in graph capture data from Pataky paper (Figure 7c - parametric results)
	# Pataky TC, Vanrenterghem J, Robinson MA. Zero- vs. one-dimensional, parametric vs
	# 	non-parametric, and confidence interval vs. hypothesis testing procedures in 
	#	one-dimensional biomechanical trajectory analysis. J Biomech 2015;48:1277-85. 
	#	doi:10.1016/j.jbiomech.2015.02.051.
	
	tmpDat = readLines('patakyT.txt');

	times = c();
	t_values = c();
	j = 1;
	for (i in 2:length(tmpDat)) {
		tmp = unlist(strsplit(tmpDat[i],"\t"));
		
		times[j] = as.numeric(tmp[1]);
		t_values[j] = as.numeric(tmp[2]);
		j = j+1;
	}
	
	tOut = data.frame(times, t_values);
	names(tOut)[1] = "Time";
	names(tOut)[2] = "tval";
	
	return (tOut);
	
}

# Last edit: Fall 2016
readPatakyData <- function() {
	# Read in data used for Pataky et al, 2015 paper
	# Pataky TC, Vanrenterghem J, Robinson MA. Zero- vs. one-dimensional, parametric vs. 
	#	non-parametric, and confidence interval vs. hypothesis testing procedures in 
	#	one-dimensional biomechanical trajectory analysis. J Biomech 2015;48:1277-85. 
	#	doi:10.1016/j.jbiomech.2015.02.051. 
	# These data are available publically as supplemental materials 
	# Besier TF, Fredericson M, Gold GE, Beaupré GS, Delp SL. Knee Muscle Forces during 
	#	Walking and Running in Patellofemoral Pain Patients and Pain-Free Controls. 
	#	J Biomech 2009;42:898-905.
	# in the file KinematicsEMGmomentarmForces-walking.csv
	 
	# This next line should be adjusted based on where YOU put these data.
	kinEMGDat = readLines('KinematicsEMGmomentarmForces-walking.csv');
	
	Subject = c();
	Group = c();
	Gender = c();
	Time = c();
	
	HipFlexion = c();
	KneeFlexion = c();
	AnkleDorsiFlexion = c();
	
	SemiMemEMG = c()
	BicepsFemEMG = c();
	RectusFemEMG = c();
	VastusMedEMG = c();
	VastusLatEMG = c();
	MedGastrocEMG = c();
	LatGastrocEMG = c();
	
	SemiMemMA = c();
	SemiTenMA = c();
	BicepsFemLH_MA = c();
	BicepsFemSH_MA = c();
	RectusFemMA = c();
	VastusMedMA = c();
	VastusIntMA = c();
	VastusLatMA = c();
	MedGastrocMA = c();
	LatGastrocMA = c();
	
	SemiMemF = c();
	SemiTenF = c();
	BicepsFemLH_F = c();
	BicepsFemSH_F = c();
	RectusFemF = c();
	VastusMedF = c();
	VastusIntF = c();
	VastusLatF = c();
	MedGastrocF = c();
	LatGastrocF = c();
	
	
	j = 1;
	for (i in 1:length(kinEMGDat)) {
		ind = i%%106;
		if (ind == 1) {
			tmp = unlist(strsplit(kinEMGDat[i],","));
			subject = tmp[2];
			group = tmp[3];
			if (group != "Control") {
				group = "Pain";	 # replace "Patellofemoral Pain"
			}
			gender = tmp[4];
			t = 1;
		}
		if (ind > 4 & ind <105) {
			tmp = unlist(strsplit(kinEMGDat[i],","));
			
			Subject[j] = subject;
			Group[j] = group;
			Gender[j] = gender;
			Time[j] = t;
			
			HipFlexion[j] = as.numeric(tmp[2])
			KneeFlexion[j] = as.numeric(tmp[3])
			AnkleDorsiFlexion[j] = as.numeric(tmp[4])
			
			SemiMemEMG[j] = as.numeric(tmp[6])
			BicepsFemEMG[j] = as.numeric(tmp[7])
			RectusFemEMG[j] = as.numeric(tmp[8])
			VastusMedEMG[j] = as.numeric(tmp[9])
			VastusLatEMG[j] = as.numeric(tmp[10])
			MedGastrocEMG[j] = as.numeric(tmp[11])
			LatGastrocEMG[j] = as.numeric(tmp[12])
			
			SemiMemMA[j] = as.numeric(tmp[14])
			SemiTenMA[j] = as.numeric(tmp[15])
			BicepsFemLH_MA[j] = as.numeric(tmp[16])
			BicepsFemSH_MA[j] = as.numeric(tmp[17])
			RectusFemMA[j] = as.numeric(tmp[18])
			VastusMedMA[j] = as.numeric(tmp[19])
			VastusIntMA[j] = as.numeric(tmp[20])
			VastusLatMA[j] = as.numeric(tmp[21])
			MedGastrocMA[j] = as.numeric(tmp[22])
			LatGastrocMA[j] = as.numeric(tmp[23])
			
			SemiMemF[j] = as.numeric(tmp[25])
			SemiTenF[j] = as.numeric(tmp[26])
			BicepsFemLH_F[j] = as.numeric(tmp[27])
			BicepsFemSH_F[j] = as.numeric(tmp[28])
			RectusFemF[j] = as.numeric(tmp[29])
			VastusMedF[j] = as.numeric(tmp[30])
			VastusIntF[j] = as.numeric(tmp[31])
			VastusLatF[j] = as.numeric(tmp[32])
			MedGastrocF[j] = as.numeric(tmp[33])
			LatGastrocF[j] = as.numeric(tmp[34])
			
			j = j+1;
			t = t+1;
		}
	}


	
	dataOut = data.frame( Subject, Group, Gender, Time, 
						HipFlexion, KneeFlexion, AnkleDorsiFlexion, 
						SemiMemEMG, BicepsFemEMG, RectusFemEMG, VastusMedEMG, VastusLatEMG, MedGastrocEMG, LatGastrocEMG, 
						SemiMemMA, SemiTenMA, BicepsFemLH_MA, BicepsFemSH_MA, RectusFemMA, VastusMedMA, VastusIntMA, VastusLatMA, MedGastrocMA, LatGastrocMA, 
						SemiMemF, SemiTenF, BicepsFemLH_F, BicepsFemSH_F, RectusFemF, VastusMedF, VastusIntF, VastusLatF, MedGastrocF, LatGastrocF);

	dataOut = dataOut[order(Time),];
	return (dataOut);
}


########## GRAPHICAL ABSTRACT #####################33
# Need to update smoothing code (May 21, 2019)
# Last edit: May 22, 2019
graphicalAbstract <- function(inputData,RFTtvalues) {
	# Outputs
	#
	#	Eight panel (actually 7 + 1 blank) graph showing RFT results in comparison to 
	#		LAAST results.  
	#
	# Inputs
	#
	#	inputData = data frame (Pataky Data) which has at least the following columns
	#		'Time','MedGastrocF', and 'Group'
	#	RFTtvalues = data frame containing t-values from rft1d (python) analysis of 
	#		same data.  This is over plotted in the t-value plot for comparison to LAAST.
	#		The columns 'Time' and 'tval' must exist for this to work.
	#		(patakyT)
	
	#library(fields);	# needed for qsreg - removed in preference of findSpan()
	library(ggplot2);
	library(egg);
	
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
	
	# 4 dof (5 curves) in each group, could be more if large bins
	# 16 and 26 
	
		s1 = findSpan(x1,y1);
		s2 = findSpan(x2,y2);
		span = min(c(s1,s2));
		print(paste('span = ',span));
	
	# Set up graphical parameters
	tRa = c(1,6.5);	# Ranges of t, dof, and p plots
	dfRa = c(20,42);
	pRa = c(-18,0);

	showX = c(T,T,T,T);	# which plots to show scales and legends on
	showY = c(T,T,T,T);
	showL = c(F,F,F,T);
	lPosMean = c(0.2,0.7);	# legend positions
	#lPosP = c(0,0);		
	lPosP = c(0.3,0.25);
	#lPosP = c(0.1,0.35);	
	showPts = F;		# show points on mean plot or no
	lHoriz = F;	# horizontal p legend
	
	gParams = list(0,100,20,"% Cycle","Musle Force (N)",span,"MedGast",tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz)
	
	# Do pooled t-test
	p1 = laastVsRFT(inputData,'Time','MedGastrocF',c('Pain','Control'),2,'pp',gParams,RFTinfo,RFTtvalues);
	#p1 = laast(inputData,'Time','MedGastrocF',c('Pain','Control'),1,'pp',gParams);
	
	p1a = p1[[1]]#+annotate("text",x=-3,y=1150,label="(a)");
	p2a = p1[[2]]#+annotate("text",x=-3,y=40,label="(c)");
	p3a = p1[[3]]+#annotate("text",x=-3,y=5.75,label="(b)")+
			annotate("text",x=20,y=2.5,label="LAAST",color="#10BE44")+
			annotate("text",x=10,y=4.5,label="RFT",color="blue")+
			annotate("text",x=35,y=6.5,label="Pooled variance t-tests",color="black")
	p4a = p1[[4]]#+annotate("text",x=-3,y=-1,label="(d)");

	showX = c(F,T,T,T);
	showY = c(F,T,T,F);
	showL = c(T,F,F,T);
	lPosMean = c(0,0);
	lPosP = c(0.3,0.25);
	#lPosP = c(0,0);
	showPts = F;
	lHoriz = F;	# horizontal p legend
	
	gParams = list(0,100,20,"% Cycle","Knee Flexion (deg)",span,"MedGast",tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz)
	
	# Do Welch t-test
	p2 = laastVsRFT(inputData,'Time','MedGastrocF',c('Pain','Control'),1,'p',gParams,RFTinfo,RFTtvalues);
	#p2 = laast(inputData,'Time','MedGastrocF',c('Pain','Control'),1,'p',gParams);
	p1b = p2[[1]]#+annotate("text",x=-3,y=1150,label="");
	p2b = p2[[2]]#+annotate("text",x=-3,y=40,label="(f)");
	p3b = p2[[3]]+#annotate("text",x=-3,y=5.75,label="(e)")+
			annotate("text",x=10,y=5.5,label="LAAST",color="#10BE44")+
			annotate("text",x=20,y=2.5,label="RFT",color="blue")+
			annotate("text",x=30,y=6.5,label="Welch's t-tests",color="black");;
	p4b = p2[[4]]#+annotate("text",x=-3,y=-1,label="(g)");
	
	# Make a blank plot
	df1 = data.frame();
	pb = ggplot(df1) + geom_point() + 
		xlim(0,5)+ ylim(0,5) + 
		annotate("text",x=2.2,y=3.5,label="LOESS Alpha-Adjusted",color="black",size=6)+
		annotate("text",x=2.2,y=2.5,label="Serial T-tests",color="black",size=6)+
		theme(line=element_blank(),text=element_blank(),panel.background = element_rect(fill = "white"))

	g1 <- ggarrange(p1a,p3a,p4a, pb,p3b,p4b, nrow=2)	# use egg rather than multiplot or other means for multiple graph
	print(g1);
	
	fN = "graphicalAbstract.png";
	ggsave(filename=fN,plot=g1,width=15,height=6,dpi=300);  # units="in",

}


##### FIGURE 1 ########
# Last edit: May 22, 2019
buildFigure1 <- function(testDF2) {
	# This function was used to create Figure 1 and is mainly a plotting function.  
	# Most of the statistical work is done in runTestCase()
	#
	# testDF2 is the output of buildData1a()
	#
	
	library(egg);
	library(ggplot2);
	
	p1 = runTestCase(testDF2,2,'Sine',isLong=FALSE);	# Sin wave
	p2 = runTestCase(testDF2,2,'Sine',isLong=TRUE);

	# How to extract graphs
	p1a = p1[[1]]+annotate("text",x=30,y=4,label="(a)");
	p2a = p1[[2]]+annotate("text",x=30,y=50,label="(b)");
	p3a = p1[[3]]+annotate("text",x=30,y=5,label="(c)");
	p4a = p1[[4]]+annotate("text",x=30,y=0,label="(d)");
	p5a = p1[[5]];
	p6a = p1[[6]]
								
	# How to extract graphs
	p1b = p2[[1]]+annotate("text",x=30,y=4,label="(e)");
	p2b = p2[[2]]+annotate("text",x=30,y=50,label="(f)");
	p3b = p2[[3]]+annotate("text",x=30,y=5,label="(g)");
	p4b = p2[[4]]+annotate("text",x=30,y=0,label="(h)");
	p5b = p2[[5]];
	p6b = p2[[6]];

	
	g1 <- ggarrange(p1a, p1b,
					p2a, p2b,
					p3a, p3b,
					p4a, p4b,
					nrow =4);	# from egg - better for multiplots
	print(g1);
	
	fN = "Figure1.png";
	ggsave(filename=fN,plot=g1,width=16,height=10,dpi=400);  # units="in",
		
}


# Called from buildFigure1() to produce Figure 1 in paper
# Last edit: May 22, 2019
# NB! Run testRFT.py first on data generated by buildData1a().  Then fill in results below.
runTestCase <- function(testDF2,span,test,isLong) {
	# testDF2 = the data frame upon which to run this test... should be of format df(Group,X,Y)
	#		Group should be either "Group1" or "Group2"
	# span = initial smoothing span.  If span < 1, then use input value and do not search for
	#		a better span to fit the data.
	# test = "Gaussian","GrowthVelocity","SineWave"
	#		The type of test determines which RFT dof and critical t-value is displayed in the graphs
	# isLong = boolean.  If TRUE, then presume longitudinal data, and sample (nsamp) without replacement.
	#
	
	if (test == 'Gaussian') {
		# RFTinfo = c(t-crit, dof)
		RFTinfo = c(2.89,38);	# Results from RFT with Gaussian (using rft1d Python code)
		tit = 'TestGaussian';
		xMax = 100;
		dX = 20;
		tRa = c(-10,0);	# Ranges of t, dof, and p plots
		dfRa = c(0,50);
		pRa = c(-20,0);
		binSize = 1;
		nsamp = 90;
		lPosMean = c(0.8,0.8);	# legend positions
	} else if (test == 'GrowthVelocity') {
		# RFTinfo = c(t-crit, dof)
		RFTinfo = c(2.987,38);	# Results from RFT with Growth Velocity simulation (using rft1d Python code)
		tit = 'TestGrowthVelocity';
		xMax = 100;
		dX = 20;
		tRa = c(-15,5);	# Ranges of t, dof, and p plots
		dfRa = c(0,50);
		pRa = c(-45,5);
		binSize = 1;
		nsamp = 90;
		lPosMean = c(0.5,0.8);	# legend positions		
	} else {					
		# RFTinfo = c(t-crit, dof)
		RFTinfo = c(2.219,38); # Results from RFT with Sine Wave (using testRFT.py Python code)
		tit = 'TestSineWave';
		xMax = 360;
		dX = 60;
		tRa = c(-5,5);	# Ranges of t, dof, and p plots
		dfRa = c(0,50);
		pRa = c(-15,0);
		binSize = 1;
		nsamp = 300;
		lPosMean = c(0.8,0.8);	# legend positions		
	}
	
	RFTtvalues = c();	# Empty vector to keep from plotting RFT t-values

	showX = c(F,F,F,T);	# which plots to show scales and legends on
	showY = c(T,T,T,T);
	showL = c(T,F,F,T);
	
	lPosP = c(0.5,0.15);			
	showPts = T;		# show points on mean plot or no
	lHoriz = T;	# horizontal p legend
	
	if (isLong) {
		ind1 = which(testDF2$Group == 'Group1');
		ind2 = which(testDF2$Group == 'Group2');

		x1 = testDF2$X[ind1];
		x2 = testDF2$X[ind2];
		y1 = testDF2$Y[ind1];
		y2 = testDF2$Y[ind2];
		
		l = length(x1)/xMax;	# number of "subjects"
		
		x1out = c();
		y1out = c();
		x2out = c();
		y2out = c();
		
		for (i in c(1:l)) {
			imin = (i-1)*xMax;
			imax = i*xMax;
			
			s1 = seq(imin,imax,1);
			ind1 = sample(s1,nsamp, replace=F);
			ind2 = sample(s1,nsamp, replace=F);
			
			X1 = x1[ind1];
			X2 = x2[ind2];
			Y1 = y1[ind1];
			Y2 = y2[ind2];
		
			x1out = c(x1out,X1);
			x2out = c(x2out,X2);
			y1out = c(y1out,Y1);
			y2out = c(y2out,Y2);
			
		}
		
		g1 = rep('Group1',length(x1out));
		g2 = rep('Group2',length(x2out));
		
		X = c(x1out,x2out);
		Y = c(y1out,y2out);
		Group = c(g1,g2);

		if (span > 1) {
			s1 = findSpan(x1out,y1out);
			s2 = findSpan(x2out,y2out);
			span = min(c(s1,s2));
			print(paste('span = ',span));
		}
		df3 = data.frame(X,Y,Group);
		gParams = list(0,xMax,dX,"Time Point","Amplitude",span,tit,
				tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz);
		p1 = laast(df3,'X','Y',c('Group1','Group2'),binSize,'p',gParams);
		
	} else {

		if (span > 1) {	# Automatically determine span
			print('test1');
			ind1 = which(testDF2$Group == 'Group1');
			print('test2');
			ind2 = which(testDF2$Group == 'Group2');
			x1 = testDF2$X[ind1];
			x2 = testDF2$X[ind2];
			y1 = testDF2$Y[ind1];
			y2 = testDF2$Y[ind2];		

			s1 = findSpan(x1,y1);
			s2 = findSpan(x2,y2);
			span = (s1+s2)/2;
			print(paste('span = ',span));

		}
		gParams = list(0,xMax,dX,"Time Point","Amplitude",span,tit,
				tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz)
		
		p1 = laastVsRFT(testDF2,'X','Y',c('Group1','Group2'),binSize,'p',gParams,RFTinfo,RFTtvalues);
	}
	
	return (p1);
}



##### FIGURE 2 ########

# Last edit: May 30, 2019
buildFigure2 <- function(inputData,RFTtvalues) {
	# This function was used to create Figure 2
	# Outputs
	#
	#	Eight panel (actually 7 + 1 blank) graph showing RFT results in comparison to 
	#		LAAST results.  
	#
	# Inputs
	#
	#	inputData = data frame (Pataky Data) which has at least the following columns
	#		'Time','MedGastrocF', and 'Group'.  These come from readPatakyData().
	#	RFTtvalues = data frame containing t-values from rft1d (python) analysis of 
	#		same data.  This is over plotted in the t-value plot for comparison to LAAST.
	#		The columns 'Time' and 'tval' must exist for this to work.
	#		(patakyT)
	#		In the paper, these data were scanned directly from Figure 7c in *********
	#		and then read into R using function -> readPatakyT()
	
	#library(fields);	# needed for qsreg - no longer used.
	
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
	
	p1a = p1[[1]]+annotate("text",x=-3,y=1150,label="(a)");
	p2a = p1[[2]]+annotate("text",x=-3,y=40,label="(c)");
	p3a = p1[[3]]+annotate("text",x=-3,y=5.75,label="(b)")+
			annotate("text",x=20,y=2.5,label="LAAST",color="#10BE44")+
			annotate("text",x=10,y=4.5,label="RFT",color="blue");
	p4a = p1[[4]]+annotate("text",x=-3,y=-1,label="(d)");

	showX = c(F,F,F,T);
	showY = c(F,F,F,F);
	showL = c(T,F,F,T);
	lPosMean = c(0,0);
	lPosP = c(0.5,0.13);
	#lPosP = c(0,0);
	showPts = F;
	lHoriz = T;	# horizontal p legend
	
	# Graphical parameters for LAAST
	gParams = list(0,100,20,"% Cycle","Knee Flexion (deg)",span,"MedGast",tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz)
	
	# Do Welch t-test
	p2 = laastVsRFT(inputData,'Time','MedGastrocF',c('Pain','Control'),1,'p',gParams,RFTinfo,RFTtvalues);
	
	p1b = p2[[1]]+annotate("text",x=-3,y=1150,label="");
	p2b = p2[[2]]+annotate("text",x=-3,y=40,label="(f)");
	p3b = p2[[3]]+annotate("text",x=-3,y=5.75,label="(e)")+
			annotate("text",x=10,y=5.5,label="LAAST",color="#10BE44")+
			annotate("text",x=20,y=2.5,label="RFT",color="blue");;
	p4b = p2[[4]]+annotate("text",x=-3,y=-1,label="(g)");
	
	# Make a blank plot
	df1 = data.frame();
	pb = ggplot(df1) + geom_point() + 
		xlim(0,5)+ ylim(0,5) + 
		theme(line=element_blank(),text=element_blank(),panel.background = element_rect(fill = "white"))

	#g1 <- grid.arrange(p1a,pb, p3a,p3b, p2a,p2b, p4a,p4b, nrow =4);
	g1 <- ggarrange(p1a,pb, p3a,p3b, p2a,p2b, p4a,p4b, nrow =4);	# from egg - better for 
	print(g1);
	
	fN = "Figure2.png";
	ggsave(filename=fN,plot=g1,width=8,height=8,dpi=400);  # units="in",

}


######### FIGURE B.2 ####################
alphaPatterns <- function() {
	# This function was used to create Figure B1
	rho = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975);
	N = c(1:50)*20;
	
	alphas = c();
	rhoArr = c();
	nArr = c();
	
	cnt = 1;
	for (r in rho) {
		for (n in N) {
			Neff = n*(1-r*r)+r*r;
			alphas[cnt] = 0.05/Neff;
			rhoArr[cnt] = r;
			nArr[cnt] = n;
			
			cnt = cnt+1;
		}
	}
	
	Rho = as.factor(rhoArr);
	plotArr = data.frame(nArr,alphas,Rho);
	
	ggplot(plotArr,aes(x=nArr, y=alphas, color=Rho))+geom_line(size=0.75)+
		xlab("Number of Comparisons")+ylab("Adjusted Alpha Level");
			
	
}


######### FIGURE B.1 & B.3 ####################
# Updated May 21, 2019
runBinDependenceBumpy <- function() {
	# This function was used to create Figure B3
	
	testDF <<- buildTestCase3a(360);

	#### FIGURE B.3  ####
	# Set up graphical parameters for Supplementary Figure B.1
	tRa = c(-6,6);	# Ranges of t, dof, and p plots
	dfRa = c(4,10);
	pRa = c(-10,0);

	showX = c(T,T,T,T);	# which plots to show scales and legends on
	showY = c(T,T,T,T);
	showL = c(T,T,T,T);
	lPosMean = c(0.6,0.7);	# legend positions
	lPosP = c(0,0);			
	showPts = T;		# show points on mean plot or no
	span = 0.33;		# Change this as needed to model sine waves	
	
	# Reinitialize these globals
	alphaNarr <<- c();
	alphaHarr <<- c();
	alphaHBarr <<- c();
	alphaBarr <<- c();
	
	# Graphical parameters for LAAST
	gParams = list(0,360,30,"Degrees","Amplitude",span,"TestSineWaveBins9",
			tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts)
	
	# Note only one of these plots is selected for the appendix.  They exhibit some variability
	# dependent on the random sample generated.
	print('test1');
	p1 = illustrateBinDependence(testDF,'X','Y',c('Group1','Group2'),9,'p',gParams);	#9
	gParams[7] = "TestSineWaveBins6"
	print('test2');
	p2 = illustrateBinDependence(testDF,'X','Y',c('Group1','Group2'),6,'p',gParams);	#6
	gParams[7] = "TestSineWaveBins3"
	print('test3');
	p3 = illustrateBinDependence(testDF,'X','Y',c('Group1','Group2'),3,'p',gParams);	#3
	gParams[7] = "TestSineWaveBins1"
	print('test4');
	p4 = illustrateBinDependence(testDF,'X','Y',c('Group1','Group2'),1,'p',gParams);	#1
	
	lab1 = rep("LAAST",16);
	lab2 = rep("Holm",16);
	lab3 = rep("Hochberg",16);
	lab4 = rep("Bonferroni",16);
	Correction = c(lab1,lab2,lab3,lab4);
	alphaBarr = 0.05/samplesArr;
	alphas = c(alphaNarr,alphaHarr,alphaHBarr,alphaBarr);
	N = c(samplesArr,samplesArr,samplesArr,samplesArr);
	
	alphaSum = data.frame(N,alphas,Correction);
	
	
	#### FIGURE B.1 ####
	p1 = ggplot(alphaSum,aes(x=N,y=alphas,color=Correction))+geom_line(size=1)+geom_point()+
		xlab("Number of Comparisons")+ylab("Adjusted Alpha");
	print(p1);
	fN = "FigureB.1.png";
	
	ggsave(filename=fN,plot=p1,width=8,height=6,dpi=300);

	return (alphaSum);
}

alphaNarr <<- c();
alphaHarr <<- c();
alphaHBarr <<- c();
samplesArr <<- c();


# Updated May 21, 2019
illustrateBinDependence <- function(allData,varX,varY,groups,binSize0,type,gParams) {
	# Start with a small binSize0 as this gets multiplied by 2, 5, and 10 
	# for comparison
	# 
	# type = 'p' which is usual, 'pp', 't', 'f', and 'tp'.  'pp' produces equal variance t-test.
	#
	# Graphical parameters for plotting
	# gParams = c(xmin,xmax,deltaX,xtit,ytit,span0,fn);
	#	where xmin,xmax are the min and max x-values that should be displayed
	#	deltaX is the x-axis tic spacing
	#	xtit is the x-axis title
	# 	ytit is the y-axis title
	# 	span0 is the span needed in the density function in nPerSamp
	#		default is 0.75 - good for data with spacing, for fine data, use smaller
	# 		values like 0.1 or 0.2.  Use 0.1 for Pataky validation
	# 	fn is the fileleaf of the file to be saved with the graphics

	library(ggplot2);
	library(gridExtra);	# for multipanel plotting with ggplot2!
	
	#plot2Screen = TRUE;
	plot2Screen = FALSE;
	
	a1Cnt <<- 1;	# count of alpha-level values from Niiler method
	a2Cnt <<- 1;	# count of alpha-level values from Holm-Bonferroni method
	
	alphaN <<- c();
	alphaH <<- c();
	alphaHB <<- c();
	
	#(dataS, varX, varY, group1, group2, gParams,binSize0)
	p1 = meanPlot(allData, varX, varY, groups[1],groups[2],gParams,binSize0);
	
	if (type == 'p' | type == 't') {
		type1 = 't';
	} else if (type == 'f') {
		type1 = 'f';
	} else {
		type1 = 'tp'
	}
	
	samples = c(binSize0,2*binSize0,4*binSize0,5*binSize0);
	
	print(type1);
	#smoothAndFit <- function(dataS, varX,varY,group1,group2,binSize0,type,gParams) {
	paramsA = smoothAndFit(allData,varX,varY,groups[1],groups[2],binSize0,type1,gParams);	# smoothed t-test
	paramsB = smoothAndFit(allData,varX,varY,groups[1],groups[2],2*binSize0,type1,gParams);	# smoothed t-test
	paramsC = smoothAndFit(allData,varX,varY,groups[1],groups[2],4*binSize0,type1,gParams);	# smoothed t-test
	paramsD = smoothAndFit(allData,varX,varY,groups[1],groups[2],5*binSize0,type1,gParams);	# smoothed t-test
	

	dfA = getParams(paramsA,type);
	dfB = getParams(paramsB,type);
	dfC = getParams(paramsC,type);
	dfD = getParams(paramsD,type);

	# where dfX = data.frame(times,test,samples,sig,degF);
	dfAll = rbind(dfA, dfB);
	dfAll = rbind(dfAll, dfC);
	dfAll = rbind(dfAll, dfD);
	
	dfAll$samples = as.factor(dfAll$samples);
	
	p2 = dofPlot(dfAll,gParams);
	#p2 = powerPlot(dfAll,gParams);
	
	sig = dfAll$sig
	
	
	# for t graph only
	if (type == 'p' | type == 'pp') {	
		int = 0.05;
		yra = c(min(sig),max(sig));
	} else {
		int = 0;
		yra = c(0,max(sig));
	}
	
	xmin = unlist(gParams[1]);
	xmax = unlist(gParams[2]);
	deltaX = unlist(gParams[3]);
	xtit = unlist(gParams[4]);
	ytit = unlist(gParams[5]);
	br = seq(xmin,xmax,deltaX);
	
	
	alphaNo = signif(0.05,digits = 3);



	if (type=='p' | type == 'pp') {
		ytit = 't-values';
	} else if (type=='t' | type == 'tp') {
		ytit = 't-values';
	} else{
		ytit = 'F-values';
	}
	
	
	p3 = ggplot(dfAll)+geom_line(aes(x=times,y=sig,color=samples), size=1.25  )+
		#geom_abline(intercept=int,slope=0)+
		ylab(ytit)+
		scale_x_continuous(breaks=br,name=xtit)+
		ylim(yra[1],yra[2])

	# Now get p-values
	if (type == 'p' | type == 't' | type == 'f') {
		type1 = 'p';
		int = log(0.05);
		yra = c(min(sig),max(sig));
		if (type == 't' | type == 'f') {
			int = 0;
			yra = c(0,max(sig));
		}
	} else {
		type1 = 'pp'
		int = log(0.05);
		yra = c(min(sig),max(sig));	
	}

	
	paramsA = smoothAndFit(allData,varX,varY,groups[1],groups[2],binSize0,type,gParams);	# smoothed t-test
	paramsB = smoothAndFit(allData,varX,varY,groups[1],groups[2],2*binSize0,type,gParams);	# smoothed t-test
	paramsC = smoothAndFit(allData,varX,varY,groups[1],groups[2],4*binSize0,type,gParams);	# smoothed t-test
	paramsD = smoothAndFit(allData,varX,varY,groups[1],groups[2],5*binSize0,type,gParams);	# smoothed t-test

	dfA = getParams(paramsA,type);
	dfB = getParams(paramsB,type);
	dfC = getParams(paramsC,type);
	dfD = getParams(paramsD,type);

	# where dfX = data.frame(times,test,samples,sig,degF);
	dfAll = rbind(dfA, dfB);
	dfAll = rbind(dfAll, dfC);
	dfAll = rbind(dfAll, dfD);

	dfAll$samples = as.factor(dfAll$samples);		
	samples = 360/c(binSize0,2*binSize0,4*binSize0,5*binSize0);
		
	CIlower = dfA$CIlower;
	CIupper = dfA$CIupper;
	times = dfA$times;
	deltaM = dfA$deltaM;
	
	CI = data.frame(times,deltaM,CIlower,CIupper);

	ytit = 'log(p-value)';
		int = log(0.05);	# log(0.05);
		print(paste("LAAST alphas",a1Cnt));
		print(alphaN);
		print(paste("Holm-Bonferroni alphas",a2Cnt));
		print(alphaH);
		print(paste("Hochberg alphas",a2Cnt));
		print(alphaHB);

			alphaNarr <<- c(alphaNarr,alphaN);
			alphaHarr <<- c(alphaHarr,alphaH);
			alphaHBarr <<- c(alphaHBarr,alphaHB);
			samplesArr <<- c(samplesArr,samples);
			
		int1 = log(alphaN[1]);
		int2 = log(alphaN[2]);
		int3 = log(alphaN[3]);
		int4 = log(alphaN[4]);
		
		int1a = log(alphaHB[1]);
		int2a = log(alphaHB[2]);
		int3a = log(alphaHB[3]);
		int4a = log(alphaHB[4]);
		
		
		print(paste(int1,int2,int3,int4));
		
		yra = c(min(sig),0);
	
	p4 = ggplot(dfAll)+geom_line(aes(x=times,y=sig,color=samples), size=1.25  )+
		geom_abline(intercept=int,slope=0)+ylab(ytit)+

		geom_abline(intercept=int1,slope=0,color="#C77CFF",linetype=2)+
		geom_abline(intercept=int2,slope=0,color="#00BFC4",linetype=2)+
		geom_abline(intercept=int3,slope=0,color="#7CAE00",linetype=2)+
		geom_abline(intercept=int4,slope=0,color="#F8766D",linetype=2)+

		geom_abline(intercept=int1a,slope=0,color="#C77CFF",linetype=3)+
		geom_abline(intercept=int2a,slope=0,color="#00BFC4",linetype=3)+
		geom_abline(intercept=int3a,slope=0,color="#7CAE00",linetype=3)+
		geom_abline(intercept=int4a,slope=0,color="#F8766D",linetype=3)+
		
		scale_x_continuous(breaks=br,name=xtit);
	
	
	p5 = ggplot(CI) + geom_line(aes(x=times,y=deltaM),size=1.25,color="black") +
		geom_line(aes(x=times,y=CIlower),size=1.25,color="red") +
		geom_line(aes(x=times,y=CIupper), size=1.25, color="red") + 
		geom_abline(intercept=0,slope=0,linetype=2,colour="black")+
		ylab('95% CI')+
		scale_x_continuous(breaks=br,name=xtit);
		
	#multiplot(p1,p2,p3,p4, cols=2);
	
		
	if (plot2Screen) {
		multiplot(p1,p2,p3,p4, cols=2);
		#multiplot(p3,p4, cols=2);
	} else {
		
		fn = unlist(gParams[7]);
		fN = paste(fn,".tiff",sep="");
		print(fN);
		#tiff(filename = fN, pointsize =12, bg = "white", res = 500)
		grid.arrange(p1, p3, p2, p4, nrow=2);
		#multiplot(p1,p2,p3,p4, cols=2);
		g <- arrangeGrob(p1,p3,p2,p4,  nrow=2);
		
		ggsave(filename=fN,plot=g,width=10,height=6,dpi=300);  # units="in",
		#dev.off();
	}
	
	plotsOut <- list(p1,p2,p3,p4,p5,p6);
	
	return (plotsOut);
}

# Called from buildFigure2  May 21, 2019
# Show a comparison of methods.  This is essentially the laast() function but with some
# extras to better highlight the comparison.
laastVsRFT <- function(allData,varX,varY,groups,binSize,type,gParams,RFTinfo,RFTtvalues) {
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
	#			tRa,dfRa,pRa,showX,showY,showL,lPosMean,lPosP,showPts,lHoriz);
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
	#	RFTinfo = vector showing results of RFT analysis via rft1d (python) c(otherTval,dfRFT);
	#		otherTval = critical t-value from this analysis
	#		dfRFT = degrees of freedom used in prior analysis
	#	RFTtvalues = data frame containing t-value results of RFT analysis via rft1d (python)
	#		and having columns 'Time' and 'tval'

	library(ggplot2);
	library(gridExtra);	# for multipanel plotting with ggplot2
	library(egg);	# Also for multipanel plotting with ggplot2 
	
	# Define globals (aids in limiting number of parameters passed to lower level functions
	
	a1Cnt <<- 1;	# count of alpha-level values from LAAST method
	a2Cnt <<- 1;	# count of alpha-level values from Holm-Bonferroni method
	
	alphaN <<- c();	# reinitialize
	alphaH <<- c();
	alphaHB <<- c();
	
	# Plot the mean
	p1 = meanPlot(allData, varX, varY, groups[1],groups[2],gParams,binSize);
	
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
	
	paramsA = smoothAndFit(allData,varX,varY,groups[1],groups[2],binSize,type1,gParams);	# smoothed t-test

	dfAll = getParams(paramsA,type);

	# where dfX = data.frame(times,test,samples,sig,degF);
	dfAll$samples = as.factor(dfAll$samples);
	
	# Plot the degrees of freedom
	p2 = dofPlot(dfAll,gParams);
	p6 = powerPlot(dfAll,gParams);
	
	sig <<- dfAll$sig
	
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
	if (class(RFTtvalues) == 'data.frame') {	# Include this in plot
		p3 = ggplot(dfAll)+geom_line(aes(x=times,y=sig), color="#10BE44",size=0.75,linetype="longdash" )+
			#geom_abline(intercept=int,slope=0)+
			ylab(ytit)+
			scale_x_continuous(breaks=br,name=xtit)+
			ylim(yra[1],yra[2])+theme3 +
			geom_line(data=RFTtvalues,aes(x=Time,y=tval),color="blue",size=0.75,linetype="solid");
	} else {
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

	dfAll$samples = as.factor(dfAll$samples);		
	
	CIlower = dfAll$CIlower;	#dfA vs dfAll?
	CIupper = dfAll$CIupper;
	times = dfAll$times;
	deltaM = dfAll$deltaM;
	
	CI = data.frame(times,deltaM,CIlower,CIupper);

	int = log(0.05);	# log(0.05);
	print(paste("LAAST alphas",a1Cnt));
	print(alphaN);
	print(paste("Holm alphas",a2Cnt));
	print(alphaH);
	print(paste("Hochberg alphas",a2Cnt));
	print(alphaHB);

	# End points of horizontal lines showing adjusted alpha levels
	int1 = log(alphaN[1]);	# LAAST alpha values
	int2 = log(alphaH[1]);	# Holm-Bonferroni alpha values
	int3 = log(alphaHB[1]); # Hochberg alpha values

	# Set up labeling for alpha value lines
	# 	Note that the order here doesn't much matter since R will automatically
	#	put legends in alphabetical order unless specifically told not to.
	#   So one just needs to get the order consistent between vectors that
	# 	comprise the dataframe legDF.
	#
	# Default values if no RFT shown
	#ints = c(int1,int2,int3);
	int1vec = rep(int1,2);	# only need end points
	int2vec = rep(int2,2);
	int3vec = rep(int3,2);
	
	lab1vec = rep("LAAST",2);
	lab2vec = rep("Holm.",2);
	lab3vec = rep("Hoch.",2);
	
	#10BE44

	lts = c("dotted", "longdash","dashed");	
	cols = c("purple","red","#20C250");	
	
	intVec = c(int1vec,int2vec,int3vec);
	labVec = c(lab1vec,lab2vec,lab3vec);
	
	xMin = min(dfAll$times);
	xMax = max(dfAll$times);	
	xV = c(xMin,xMax,xMin,xMax,xMin,xMax);
	
	#ltype = c("dashed","dotted","longdash");
	
	if (RFTinfo[1] > -1) {	# Add a third line if RFT is to be shown
		# Determine alpha level corresponding to 
		#dfRFT = getDFstd(allData,groups[1],groups[2]);
		dfRFT = RFTinfo[2];
		otherTval = RFTinfo[1];
		pRFT = 2*pt(otherTval,dfRFT,lower=FALSE);
		int4 = log(pRFT);
		print("RFT alpha");
		print(pRFT);
		
		lts = c("dotted", "longdash","dashed","solid");
		cols = c("purple","red","#20C250","blue");
		
		int4vec = rep(int4,2);
		lab4vec = rep("RFT",2);
		intVec = c(int1vec,int2vec,int3vec,int4vec);
		labVec = c(lab1vec,lab2vec,lab3vec,lab4vec);
		
		xV = c(xMin,xMax,xMin,xMax,xMin,xMax,xMin,xMax);
	}
	
	
	legDF = data.frame(xV,intVec,labVec);
	names(legDF)[3] = "Alpha"; # renaming labels
	
	# New ytit for next graph
	ytit = 'log(p-value)';
		
	yra = c(min(sig),0);
	if (max(pRa)-min(pRa) != 0) {
		yra = pRa;
	}
		
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
	
	
	# p-value plots - will be ordered alphabetically Hochberg, Holm, LAAST, RFT
	# and colors/linetypes have been chosen with this in mind
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
		
	g1 <- ggarrange(p1, p3, p2, p4, nrow =2);	# from egg - better for 
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
	
	plotsOut <- list(p1,p2,p3,p4,p5,p6);
	
	return (plotsOut);
}
