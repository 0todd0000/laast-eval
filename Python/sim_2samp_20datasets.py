
import os,pathlib
import numpy as np
from spm1d import rft1d



def write_csv(fname, y0, y1):
	y   = np.vstack([y0, y1])                    # vertically stack data
	g   = [0]*y0.shape[0] + [1]*y1.shape[0]      # group indices
	a   = np.vstack( [g, y.T] ).T                # add group indices to the array
	fmt = ('%d,' + '%.3f,'*y.shape[1])[:-1]      # string format
	np.savetxt(fname, a, delimiter=',', fmt=fmt) # write file



#(0) Generate 20 datasets:
dirREPO = pathlib.Path( __file__ ).parent.parent
dirDATA = os.path.join( dirREPO, 'Data', '20datasets')
J       = 5
Q       = 101
FWHM    = 20
np.random.seed(10)
for i in range(20):
	y0  = rft1d.randn1d(J, Q, FWHM, pad=True)
	y1  = rft1d.randn1d(J, Q, FWHM, pad=True)
	write_csv( os.path.join(dirDATA, f'{i+1}.csv') , y0, y1)



