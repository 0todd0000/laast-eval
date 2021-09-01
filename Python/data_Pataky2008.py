
# This script exports the "Pataky2008-grf.csv" data file
#   The original data are available in spm1d
#
#   The main purpose of this data file was to show that
#      the LOESS procedure as implemented in LAASTv2
#      may have undesirable biomechanical consequences
#      when analyzing data with potentially important
#      high-frequency components, like impact forces
#
#   The following script analyzes these exported data:
#      ./R/scripts/ex_loess.R


import os,pathlib
import numpy as np
from matplotlib import pyplot as plt
import spm1d


# load data
dataset      = spm1d.data.uv1d.anova1.SpeedGRFcategorical()
y,A          = dataset.get_data()
y            = y[A==3]   # fast walking


# plot:
plt.close('all')
plt.figure( figsize=(8,6) )
ax = plt.axes()
ax.plot(y.T, 'b-', lw=0.5)
ax.plot(y.mean(axis=0), 'k-', lw=5)
plt.show()



# save:
dirREPO      = pathlib.Path( __file__ ).parent.parent
fnameCSV     = os.path.join( dirREPO, 'Data', 'Pataky2008-grf.csv')
np.savetxt(fnameCSV, y, delimiter=',', fmt='%.5f')




