
# This script uses a random number generator from spm1d.rft1d
#   to create the data files:
#      Gaussian_FWHM=20_Q=100.csv
#      Gaussian_FWHM=20_Q=101.csv
#
#   The main purpose of these data files was to test how
#      LAASTv2 handles data with an arbitrary number of
#      domain nodes (Q)



import os,pathlib
import numpy as np
from matplotlib import pyplot as plt
from spm1d import rft1d



# create data:  (set Q to 100 or 101 to create the two data files in the repository)
np.random.randn(0)
J,Q,W        = 8, 100, 20
yA           = rft1d.randn1d(J, Q, W, pad=True)
yB           = rft1d.randn1d(J, Q, W, pad=True)


# plot:
plt.close('all')
ax = plt.axes()
ax.plot(yA.T, 'k', lw=0.3)
ax.plot(yB.T, 'r', lw=0.3)
ax.plot(yA.mean(axis=0), 'k', lw=3)
ax.plot(yB.mean(axis=0), 'r', lw=3)
plt.show()



# save:
dirREPO      = pathlib.Path( __file__ ).parent.parent
fnameCSV     = os.path.join( dirREPO, 'Data', f'Gaussian_FWHM=20_Q={Q}.csv')
# stack data into an array with a GROUPS column:
y            = np.vstack([yA, yB])
g            = [0]*yA.shape[0] + [1]*yB.shape[0]
a            = np.vstack( [g, y.T] ).T
# write to file:
fmt          = ('%d,' + '%.3f,'*y.shape[1])[:-1]
np.savetxt(fnameCSV, a, delimiter=',', fmt=fmt)




