
import os,pathlib
import numpy as np
from matplotlib import pyplot as plt
import statsmodels.api as sm
import spm1d
lowess = sm.nonparametric.lowess


# load data:
dirREPO      = pathlib.Path( __file__ ).parent.parent
dirDATA      = os.path.join( dirREPO, 'Data')
fnameCSV0    = os.path.join( dirDATA, 'Pataky2008-grf.csv' )
fnameCSV1    = os.path.join( dirDATA, 'Pataky2008-grf-loess-mean.csv' )
y            = np.loadtxt(fnameCSV0, delimiter=',')
q0,m0        = np.loadtxt(fnameCSV1, delimiter=',')


# plot:
plt.close('all')
plt.figure( figsize=(6,4.5) )
fontname = 'Arial'
ax = plt.axes()
h0 = ax.plot(y.T, '-', lw=0.5, color='0.5')[0]
h1 = ax.plot(y.mean(axis=0), 'k-', lw=3)[0]
h2 = ax.plot(q0-1, m0, 'ro-', lw=1)[0]
labels = 'Individual subjects', 'Mean estimate (least-squares)', 'Mean estimate (LOESS, Niiler 2020)'
leg    = ax.legend( [h0,h1,h2], labels )
plt.setp(leg.get_texts(), size=10, name=fontname)
ax.set_xlabel('Time  (%)', size=14, name=fontname)
ax.set_ylabel('Vertical GRF (BW)', size=14, name=fontname)
plt.setp( ax.get_xticklabels() + ax.get_yticklabels(), size=10, name=fontname  )
plt.tight_layout()
plt.show()


# # save:
# fnamePDF     = os.path.join( dirREPO, 'Figures', 'fig_loess.pdf' )
# plt.savefig( fnamePDF )