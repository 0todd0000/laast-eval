
import numpy as np;
import rft1d;
import random;
import matplotlib.pyplot as plt;

d = 0;	# effect size, offset for power test
n = 10000; 	# Number of curves in sample distribution

# Number of curves in each sample.  Change this to match Table C1
n1 = 20;	#
n2 = 20;    #
nNodes = 360;
x = np.arange(nNodes);
y = np.zeros((n,nNodes))
x1 = np.zeros((n,nNodes));


sampArr = np.arange(n);	#
for i in range(n):
	y[i] = np.sin(np.pi*x/180)+random.gauss(0,1);
	x1[i] = x;
	

nruns = 10000;		# Number of iterations
nullRejected = 0;

print("Starting sampling");


for j in range(nruns):
	y1 = np.zeros((n1,nNodes)) #[]; #[0]*nsamp;
	y2 = np.zeros((n2,nNodes));
		
	x1a = np.zeros((n1,nNodes));
	
	s1 = np.random.choice(sampArr,(n1), replace=False);
	for i in range(n1):
		random.seed();
		rn = s1[i];
		y1[i] = y[rn];
		x1a[i] = x;
		
	x2a = np.zeros((n2,nNodes));	
	s2 = np.random.choice(sampArr,(n2), replace=False);
	for i in range(n2):
		rn = s2[i];
		if (d != 0):
			y2[i] = y[rn]+d;
		else:
			y2[i] = y[rn];
			
		x2a[i] = x;

	x1a.sort();
	x2a.sort();
	
	y1.sort();
	y2.sort();

	m1 = y1.mean(axis=0);
	m2 = y2.mean(axis=0);

	s1 = y1.std(ddof = 1, axis=0);
	s2 = y2.std(ddof = 1, axis=0);

	df = n1 + n2 - 2;
	s = np.sqrt(((n1-1)*s1*s1 + (n2-1)*s2*s2)/df);
	t = (m1-m2)/(s*np.sqrt(1.0/n1 + 1.0/n2));
	
	r1 = y1 - m1;
	r2 = y2 - m2;

	residuals = np.vstack([r1, r2]);
	FWHM = rft1d.geom.estimate_fwhm(residuals);

	alpha = 0.025; # Two tailed

	t_critical = rft1d.t.isf(alpha, df, nNodes, FWHM)

	tmax = max(t);
	tmin = min(t);

	# Two tailed
	if abs(tmin) > abs(tmax):
		tmax = tmin;

	tmax = abs(tmax);

	if tmax > abs(t_critical):
		nullRejected = nullRejected + 1;
	
print("Two Tailed Results");
percentRejected = 100.*nullRejected/nruns;
if (d != 0):
	print("********** Power is: **************");
else:
	print("******** False Positive Rate %");
	
print(percentRejected);
print("  ");



