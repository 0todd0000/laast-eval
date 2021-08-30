
import numpy as np
from matplotlib import pyplot as plt
import spm1d



#(0) Load data:
dataset      = spm1d.data.uv1d.t2.SimulatedTwoLocalMax()
yA,yB        = dataset.get_data()  #A:slow, B:fast



y            = np.vstack([yA, yB])
g            = [0]*yA.shape[0] + [1]*yB.shape[0]
a            = np.vstack( [g, y.T] ).T
fname1       = '/Users/todd/GitHub/laast-eval/Data/SimulatedTwoLocalMax.csv'
fmt          = ('%d,' + '%.3f,'*y.shape[1])[:-1]
np.savetxt(fname1, a, delimiter=',', fmt=fmt)


# #(1) Conduct test:
# alpha        = 0.05
# T2           = spm1d.stats.hotellings2(YA, YB)
# T2i          = T2.inference(0.05)
# print( T2i )
#

#(2) Plot:
plt.close('all')
ax = plt.axes()
ax.plot(yA.T, 'k', lw=0.3)
ax.plot(yB.T, 'r', lw=0.3)
ax.plot(yA.mean(axis=0), 'k', lw=3)
ax.plot(yB.mean(axis=0), 'r', lw=3)
plt.show()




# fname = '/Users/todd/GitHub/spm1d/spm1d/data/datafiles/Besier2009muscleforces.npz'
# with np.load(fname) as Z:
# 	print( list(Z.keys()) )