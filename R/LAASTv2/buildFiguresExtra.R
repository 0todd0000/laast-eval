# buildFiguresExtra.R


##### Extra figure comparison using growth velocity simulation ######
# This code is essentially what is used in Figure 1, but with minor modifications in labeling
#
# This is there to demonstrate how the smoothing algorithm is somewhat 
# robust to sharper bumps in the data.
#
# Last edit: May 22, 2019
buildFigure1Extra <- function(testDF_6GV) {
	# testDF_6GV is created by the function buildGrowthVelocityData()

	library(egg);
		
	p3 = runTestCase(testDF_6GV,2,'GrowthVelocity',isLong=FALSE);	# Gaussian bump
	p4 = runTestCase(testDF_6GV,2,'GrowthVelocity',isLong=TRUE);

	# How to extract graphs
	p1c = p3[[1]]+annotate("text",x=10,y=1.5,label="(a)");
	p2c = p3[[2]]+annotate("text",x=10,y=50,label="(b)");
	p3c = p3[[3]]+annotate("text",x=10,y=5,label="(c)");
	p4c = p3[[4]]+annotate("text",x=10,y=5,label="(d)");
	p5c = p3[[5]];
	p6c = p3[[6]];

	# How to extract graphs
	p1d = p4[[1]]+annotate("text",x=10,y=1.5,label="(e)");
	p2d = p4[[2]]+annotate("text",x=10,y=50,label="(f)");
	p3d = p4[[3]]+annotate("text",x=10,y=5,label="(g)");
	p4d = p4[[4]]+annotate("text",x=10,y=5,label="(h)");
	p5d = p4[[5]];
	p6d = p4[[6]];

	
	g1 <- ggarrange(p1c, p1d,
					p2c, p2d,
					p3c, p3d,
					p4c, p4d,
					nrow=4);	# from egg - better for 
	print(g1);
	
	fN = "Figure1bExtra.png";
	ggsave(filename=fN,plot=g1,width=10,height=10,dpi=400);  # units="in",
		

}

# This next function is not used in the paper, but is used to generate sample growth velocity
# curves for additional testing of LAAST on odd shapes and longitudinal type data similar as to 
# what may be found in:
#	M.B. Bober, T. Niiler, A.L. Duker, J.E. Murray, T. Ketterer, M.E. Harley, 
#		S. Alvi, C. Flora, C. Rustad, E.M.H.F. Bongers, L.S. Bicknell, C. Wise, 
#		A.P. Jackson, Growth in individuals with Majewski osteodysplastic primordial 
#		dwarfism type II caused by pericentrin mutations, Am. J. Med. Genet. A. 158A 
#		(2012) 2719-2725. doi:10.1002/ajmg.a.35447.
# 
buildGrowthVelocityData <- function(f) {
	# f is a scaling factor for the second group which should range from 2-6 in increments of 1
	
	x = seq(-80,20,length=100)
	hx = dnorm(x,sd=3)

	x1 = x+80
	#plot(x1,hx,type="l",lty=2)
	
	X1 = c();
	X2 = c();
	Y1 = c();
	Y2 = c();
	
	Yout = c();

	for (i in 1:20) {
	
		Y = hx+rnorm(1,0,0.1) - (4e-6)*(x1-50)^3;
		Y1 = c(Y1,Y);
		Y = f*hx+rnorm(1,0,0.1) - (4e-6)*(x1-50)^3;
		Y2 = c(Y2,Y);
		
		X1 = c(X1,x1);
		X2 = c(X2,x1);
	}
	#X1 = sample(Xvec,N1,replace=F);
	#X2 = sample(Xvec,N2,replace=F);
	l = length(X1);
	
	name1 = rep("Group1",l);
	name2 = rep("Group2",l);
	
	# keep from going below zero
	y1a = abs(min(Y1));
	Y1 = Y1 + y1a;
	y2a = abs(min(Y2));
	Y2 = Y2 + y2a;
	
	X = c(X1,X2);
	Y = c(Y1,Y2);
	Group = c(name1, name2);
	
	testDF = data.frame(X,Y,Group);
	
	testDFa = data.frame(X1,Y1,name1);
	names(testDFa) = c("X","Y","Group");
	fn = paste("dataGaussGrowthVelocityTest",f,"a.txt",sep="");
	write.table(testDFa,fn,sep="\t",quote=FALSE,row.names=FALSE);

	testDFb = data.frame(X2,Y2,name2);
	names(testDFb) = c("X","Y","Group");
	fn = paste("dataGaussGrowthVelocityTest",f,"b.txt",sep="");
	write.table(testDFb,fn,sep="\t",quote=FALSE,row.names=FALSE);
	
	
	return (testDF);
	
}
