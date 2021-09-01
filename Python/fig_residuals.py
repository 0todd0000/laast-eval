

import os,pathlib
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import patches
from spm1d import rft1d




#(0) Load experimental data and compute residuals:
dirREPO      = pathlib.Path( __file__ ).parent.parent
dirDATA      = os.path.join( dirREPO, 'Data')
fnameCSV0    = os.path.join( dirDATA, 'Schwartz2008.csv' )
fnameCSV1    = os.path.join( dirDATA, 'Kautz1991.csv' )
y0           = np.loadtxt(fnameCSV0, delimiter=',')
y1           = np.loadtxt(fnameCSV1, delimiter=',')
r0,r1        = y0 - y0.mean(axis=0) , y1 - y1.mean(axis=0)  # residuals
q0,q1        = np.linspace(0, 100, y0.shape[1]), np.linspace(0, 100, y1.shape[1])  # domain points



#(1) Generate Gaussian noise
np.random.seed(20)
J0       = 200   # sample size used to visually represent population
J        = 8     # sample size
Q        = 101   # number of domain points
FWHM     = 20    # smoothness
y2       = rft1d.randn1d(J, Q, FWHM, pad=True)    # random sample
y20      = rft1d.randn1d(200, Q, FWHM, pad=False) # visual population representation
r2       = y2 - y2.mean(axis=0)    # residuals
q2       = np.linspace(0, 100, Q)  # domain points



#(2) Read simulated data (from Niiler, 2021 - response letter):
def read_Niiler2021(fnameCSV):
	with open(fnameCSV, 'r') as f:
		lines = f.readlines()
	a         = np.array( [s.strip().split(',')[1:]  for s in lines[1:]], dtype=float )
	x,y,subj  = a.T
	subj      = np.asarray(subj, dtype=int)
	usubj     = np.unique(subj)
	q         = np.array([x[subj==u]  for u in usubj])[0]   # domain nodes
	y         = np.array([y[subj==u]  for u in usubj])      # dependent variable
	return q, y
fnameCSV  = os.path.join( dirDATA, 'Niiler2021-Fig2Dist.csv' )
q3,y3     = read_Niiler2021(fnameCSV)
r3        = y3 - y3.mean(axis=0)
# create representations of these residuals
#    Both (1) the actual residuals, and (2) random representations
#    of those residuals are plotted in the figure this script genreates
n         = 200   # sample size for visual representation of the population
# horizontal lines with Gaussian variance:
n30       = (np.random.randn(n)* np.array( [[1,1]]*n ).T).T
# straight lines with Gaussian variance:
n31       = np.vstack([np.random.randn(n), np.random.randn(n)]).T
n3        = np.vstack( [n30,n31 ])   # approximate represenation of the residuals used in Niiler's (2021) response letter




#(2) Plot:
plt.close('all')
fig = plt.figure( figsize=(10,7.5) )

axx,axy = np.linspace(0.048, 0.77, 4), np.linspace(0.70, 0.01, 4)
axw,axh = 0.22, 0.22

AX      = np.array([[plt.axes([xx,yy,axw,axh])  for yy in axy]  for xx in axx]).T


plt.get_current_fig_manager().window.move(0, 0)
fontname = 'Arial'
fsize_big = 14
### plot Schwartz data:
joint = 3


# plot Schwartz data:
ax0,ax1,ax2,ax3 = AX[:,0]
tx0 = ax0.text(0.5, 0.5, 'Mean Unknown', transform=ax0.transAxes)
tx1 = ax1.text(0.5, 0.5, 'Variance Unknown', transform=ax1.transAxes)
h0  = ax2.plot(q0, y0.T, '0.5', lw=0.2)[0]
h1  = ax2.plot(q0, y0.mean(axis=0), 'k', lw=3)[0]
hr  = ax3.plot(q0, r0.T, 'k', lw=0.2)
ax3.axhline(0, color='k', lw=3)
plt.setp([tx0,tx1], size=fsize_big, name=fontname, ha='center', va='center')
# leg = ax2.legend([h0,h1], ['Observation', 'Sample mean'], loc='lower left', bbox_to_anchor=(0.7,0))
# plt.setp(leg.get_texts(), size=10, name=fontname)
# plt.setp(hr[5], color='r', lw=3, zorder=10)
# plt.setp(hr[28], color='r', lw=3, zorder=10)
hr0,hr1   = hr[31], hr[15]
plt.setp(hr0, color='r', lw=3, zorder=10)
# plt.setp(hr[14], color='c', lw=3, zorder=10)
plt.setp(hr1, color='c', lw=3, zorder=10)



# plot Kautz data:
ax0,ax1,ax2,ax3 = AX[:,1]
tx0 = ax0.text(0.5, 0.5, 'Mean Unknown', transform=ax0.transAxes)
tx1 = ax1.text(0.5, 0.5, 'Variance Unknown', transform=ax1.transAxes)
h0  = ax2.plot(q1, y1.T, '0.5', lw=0.5)[0]
h1  = ax2.plot(q1, y1.mean(axis=0), 'k', lw=3)[0]
hr  = ax3.plot(q1, r1.T, '0.5', lw=0.5)
ax3.axhline(0, color='k', lw=3)
plt.setp([tx0,tx1], size=fsize_big, name=fontname, ha='center', va='center')
leg0= ax2.legend([h0,h1], ['Observation', 'Sample mean'], loc='lower left', bbox_to_anchor=(-0.47,0))
leg1= ax3.legend([hr0,hr1], ['Residual for discussion', 'Residual for discussion'], loc='lower left', bbox_to_anchor=(-0.47,0))
plt.setp(leg0.get_texts() + leg1.get_texts(), size=10, name=fontname)


plt.setp(hr[2], color='r', lw=3, zorder=10)
# plt.setp(hr[9], color='c', lw=2)
# plt.setp(hr[10], color='c', lw=2)
plt.setp(hr[18], color='c', lw=3, zorder=10)

# plot RFT simulation:
ax0,ax1,ax2,ax3 = AX[:,2]
ax0.axhline(0, color='k', lw=3)
ax0.set_ylim(-3, 3)
ax1.plot( q2, y20.T, 'k', lw=0.2)
ax2.plot( q2, y2.T, '0.5', lw=0.5)
ax2.plot( q2, y2.mean(axis=0), 'k', lw=3 )
ax2.axhline(0, color='k', ls=':')
ax3.plot(q2, r2.T, '0.5', lw=0.5)
ax3.plot(q2, r2.mean(axis=0), 'k', lw=3)

plt.setp([ax2,ax3], ylim=(-3,3))


# plot LAAST simulation:
ax0,ax1,ax2,ax3 = AX[:,3]
ax0.plot( q3, y3.mean(axis=0), 'r', lw=3 )
ax0.axhline(0, color='k', ls=':')
ax0.set_ylim(-3.5, 3.5)
ax1.plot( n3.T, 'r', lw=0.2)
ax2.plot( q3, y3.T, color=(1,0.7,0.7), lw=0.2)
ax2.plot( q3, y3.mean(axis=0), 'k', lw=3 )
ax2.axhline(0, color='k', ls=':')
ax3.plot(q3, r3.T, 'r', lw=0.2)
ax3.plot(q3, r3.mean(axis=0), 'k', lw=3)


# column labels:
labels0 = 'Schwartz et al. (2008)', 'Kautz et al. (1991)', 'Liebl (2021)', 'Niiler(2021)'
labels1 = 'Hip flex/ext (+/-)', 'Normal pedal force'
txx0    = [ax.set_title(s)  for ax,s in zip(AX[0], labels0)]
txx1    = [ax.text(0.5, 0.9, s, transform=ax.transAxes)  for ax,s in zip(AX[0], labels1)]
plt.setp(txx0, name=fontname, size=fsize_big)
plt.setp(txx1, name=fontname, size=11, ha='center')


# row labels:
labels0 = 'Population Mean', 'Population Variance', 'Dataset', '`Residuals', 'Statistical Inference'
labels1 = '', '', '(Random sample)', '(Null hypothesis model)', ''
txx0    = [ax.text(-0.20, 0.5, s, transform=ax.transAxes)  for ax,s in zip(AX[:,0], labels0)]
txx1    = [ax.text(-0.12, 0.5, s, transform=ax.transAxes)  for ax,s in zip(AX[:,0], labels1)]
plt.setp(txx0 + txx1, name=fontname, size=13, rotation=90, va='center')
plt.setp(txx1, size=11)


# multi-column labels
labels = 'SPM Experimental Data Analysis', 'Numerical Simulation'
txx    = [ax.text(1.1, 1.22, s, transform=ax.transAxes)  for ax,s in zip(AX[0,[0,2]], labels)]
plt.setp(txx, size=18, ha='center', name=fontname)
[plt.setp( ax.get_xticklabels() + ax.get_yticklabels(), size=10, name=fontname)  for ax in AX.ravel()]


plt.setp(AX, xticks=[], yticks=[])


# inter-panel arrows
def plot_arrow_connection(ax0, ax1):
	arrow = patches.ConnectionPatch([0.5,0.15], [0.5,0.95], coordsA=ax0.transAxes, coordsB=ax1.transAxes,
		color="0.2", arrowstyle="-|>",  mutation_scale=20,  linewidth=2)
	fig.patches.append(arrow)
[[plot_arrow_connection(AXX[i], AXX[i+1])  for i in range(3) ]  for AXX in AX.T]


# background patch:
rect = patches.Rectangle([0.041,0], 0.476, 1, transform=fig.transFigure, zorder=-1, fc='0.7')
fig.patches.append( rect )

# panel letter labels
txx = [ax.text(0.01, 0.92, f'({chr(i+97)})', transform=ax.transAxes)  for i,ax in enumerate(AX.ravel())]
plt.setp(txx, size=11, name=fontname)

plt.show()


# save:
fnamePDF     = os.path.join( dirREPO, 'Figures', 'fig_residuals.pdf' )
plt.savefig( fnamePDF )

