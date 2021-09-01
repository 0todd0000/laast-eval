
import os,pathlib
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import spm1d


def load_csv(fname):
	a     = np.loadtxt(fname, delimiter=',')
	g,y   = a[:,0], a[:,1:]
	y0,y1 = y[g==0], y[g==1]
	return y0,y1



#(0) Load the 20 datasets and analyze using SPM:
dirREPO     = pathlib.Path( __file__ ).parent.parent
dirDATA     = os.path.join( dirREPO, 'Data', '20datasets')
alpha       = 0.05
Y0,Y1       = [], []  # datasets
T           = []      # t-value continua
tcrit_spm   = []      # critical t values
for i in range(20):
	fnameCSV  = os.path.join(dirDATA, f'{i+1}.csv')
	y0,y1     = load_csv( fnameCSV )
	ti        = spm1d.stats.ttest2(y0, y1).inference(alpha, two_tailed=True)
	tcrit_spm.append( ti.zstar )
	Y0.append( y0 )
	Y1.append( y1 )
	T.append( ti.z )



#(1) Load the LAAST results (see "./R/sim_2samp_20datasets.R")
fnameLAAST  = os.path.join( dirREPO, 'Data', '20datasets-laast-results.csv')
a           = np.loadtxt(fnameLAAST, delimiter=',', skiprows=1)
pcrit_laast = a[:,1]
tcrit_laast = stats.t.isf( pcrit_laast, 8 )   # approximate LAAST critical threhold




#(2) Plot:
plt.close('all')
fig,AX   = plt.subplots( 5, 8, figsize=(12,8) )
AX0      = AX[:,::2]   # axes for datasets
AX1      = AX[:,1::2]  # axes for t values
fontname = 'Helvetica'
color_rft   = 'red'
color_laast = (0.3,0.3,0.7)
for i,(ax0,ax1,y0,y1,t) in enumerate( zip(AX0.ravel(), AX1.ravel(), Y0, Y1, T) ):
	# plot dataset:
	ax0.plot(y0.T, 'k', lw=0.5)
	ax0.plot(y1.T, 'c', lw=0.5)
	# plot statistical results:
	ax1.plot( t, color='0.7', lw=3, label='t value' )
	ax1.axhline( 0, color='k', ls='-', lw=0.3 )
	ax1.axhline( tcrit_spm[i],   color=color_rft,   ls='--', label='RFT' )
	ax1.axhline( tcrit_laast[i], color=color_laast, ls='--', label='LAAST' )
	# dataset label:
	ax0.text(0.5, 1.1, f'Dataset {i+1}', color='k', fontweight='bold', ha='center', size=10, transform=ax0.transAxes, name=fontname)
	# H0 rejection checks:
	tmax = t.max()
	if tmax > tcrit_spm[i]:
		ax1.text(0.5, 0.20, 'Significant', color=color_rft, fontweight='bold', ha='center', size=8, transform=ax1.transAxes, name=fontname)
	if tmax > tcrit_laast[i]:
		ax1.text(0.5, 0.05, 'Significant', color=color_laast, fontweight='bold', ha='center', size=8, transform=ax1.transAxes, name=fontname)
	if (tmax > tcrit_spm[i]) or (tmax > tcrit_laast[i]):
		ax1.set_facecolor( '0.8' )
plt.setp(AX, ylim=(-5, 5))
leg0  = AX[0,0].legend( [AX[0,0].lines[0], AX[0,0].lines[-1]], ['Group A', 'Group B'], bbox_to_anchor=(0.91,0.25), loc='upper right' )
leg1  = AX[0,1].legend( bbox_to_anchor=(0.91,0.25), loc='upper right' )
plt.setp( leg0.get_texts() + leg1.get_texts(), name=fontname, fontsize=8 )
[plt.setp(ax, xticklabels=(), yticklabels=())  for ax in AX.ravel()[1:]]
[plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), size=7, name=fontname)  for ax in AX.ravel()]
plt.tight_layout()
plt.show()



# save:
fnamePDF  = os.path.join( dirREPO, 'Figures', 'fig_2samp_20datasets.pdf' )
plt.savefig( fnamePDF )
