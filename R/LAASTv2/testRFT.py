
# Path to files to test
#path1 = "dataGaussTest2a.txt"
#path2 = "dataGaussTest2b.txt"
#path1 = "dataGaussGrowthVelocityTest4a.txt"
#path2 = "dataGaussGrowthVelocityTest4b.txt"
path1 = "dataSineTestPhase2a.txt";
path2 = "dataSineTestPhase2b.txt";

import numpy as np;
import rft1d;

data1 = np.loadtxt(path1, dtype={'names': ('X', 'Y','Group'), 'formats': ('f4', 'f4', 'S6')}, skiprows=1);
data2 = np.loadtxt(path2, dtype={'names': ('X', 'Y','Group'), 'formats': ('f4', 'f4', 'S6')}, skiprows=1);


nNodes = 100;	# Corresponds to binSize0 = 1;
n = 20;		# Number of "subjects"
n1 = 20; 
n2 = 20;

print(path1);

y1 = data1['Y'];
y2 = data2['Y'];

Y1 = np.zeros((n,nNodes))
Y2 = np.zeros((n,nNodes))
for i in range(n):
	for j in range(nNodes):
		ind = (i*n)+j;
		Y1[i][j] = y1[ind];
		Y2[i][j] = y2[ind];


print(y1);

m1 = Y1.mean(axis=0);
m2 = Y2.mean(axis=0);

r1 = Y1 - m1;
r2 = Y2 - m2;

residuals = np.vstack([r1, r2]);
FWHM = rft1d.geom.estimate_fwhm(residuals);
print('FWHM = ')
print(FWHM);


alpha = 0.05
df = n1 + n2 - 2;

print("nNodes = ",nNodes);
t_critical = rft1d.t.isf(alpha, df, nNodes, FWHM)
print("t_crit = ",t_critical);

c = 2;
k = 20/FWHM;
tstar = t_critical;
Pset = rft1d.t.p_set(c, k, tstar, df, nNodes, FWHM)
print("Pset = ",Pset);


nNodes = 50;	# Corresponds to binSize0 = 2;
print("nNodes = ",nNodes);
t_critical = rft1d.t.isf(alpha, df, nNodes, FWHM)
print("t_crit = ",t_critical);
k = 10/FWHM;
tstar = t_critical;
Pset = rft1d.t.p_set(c, k, tstar, df, nNodes, FWHM)
print("Pset = ",Pset);

nNodes = 25;	# Corresponds to binSize0 = 4;
print("nNodes = ",nNodes);
t_critical = rft1d.t.isf(alpha, df, nNodes, FWHM)
print("t_crit = ",t_critical);
k = 4/FWHM;
tstar = t_critical;
Pset = rft1d.t.p_set(c, k, tstar, df, nNodes, FWHM)
print("Pset = ",Pset);


