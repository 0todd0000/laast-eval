
# This script runs a sequence of 20 two-sample t tests
#    using simple scalar values from the standard normal
#    distribution, then saves the results in the file:
#    ./Tables/table_val0d.tex
#
# The points are:
#    - Approximately 1 in 20 tests reach significance
#      when the null hypothesis is true (i.e. there is
#      truly no difference between the groups)
#
#    - In addition to the data itself, sample means,
#      sample standard deviations, and sample t values
#      all also vary. Inference procedures must 
#      robustly handle of this variance.


import os,pathlib
import numpy as np
from scipy import stats


def tstat(yA, yB):
	# two-sample t statistic (using pooled SD)
	mA,mB  = yA.mean(), yB.mean()
	sA,sB  = yA.std(ddof=1), yB.std(ddof=1)
	nA,nB  = yA.size, yB.size
	sp     = (   (  (nA-1)*sA*sA + (nB-1)*sB*sB  )  /  ( nA+nB-2 )   )**0.5
	t      = (mA-mB) / sp / (1.0/nA + 1.0/nB)**0.5
	return t



#(0) Generate datasets:
np.random.seed(7)
J       = 5     # sample size
alpha   = 0.05  # Type I error rate
tcrit   = stats.t.isf(0.05/2, 2*J-2)   # critical t statistic (two-tailed)
results = []
for i in range(20):
	y0     = np.random.randn(J)     # random data (standard normal)
	y1     = np.random.randn(J)     # random data (standard normal)
	md     = y1.mean() - y0.mean()  # mean difference
	t      = tstat(y1, y0)          # t statistic
	fp     = abs(t) > tcrit         # false positive (two-tailed)
	results.append( [y0, y1, md, t, fp] )
### print results (to screen):
print('Dataset  Mean (B-A)    t     Significant?')
for i,(y0,y1,md,t,fp) in enumerate( results ):
	print( f'{i:^8} {md:^10.3f} {t:^8.3f} {fp}' ) # i, md, t, fp )



#(1) Write table in LaTeX format
dirREPO      = pathlib.Path( __file__ ).parent.parent
dirDATA      = os.path.join( dirREPO, 'Data')
fnameTEX     = os.path.join( dirREPO, 'Tables', 'table_val0d.tex' )
# specify LaTeX frontmatter:
frontmatter = []
frontmatter.append( '\\textwidth = 450pt' )
frontmatter.append( '\\textheight = 620pt' )
frontmatter.append( '\\oddsidemargin = 10pt' )
frontmatter.append( '\\voffset = 0pt' )
frontmatter.append( '\\topmargin = 0pt' )
frontmatter.append( '\\headheight = 0pt' )
frontmatter.append( '\\footskip = 0pt' )
# write file:
with open(fnameTEX, 'w') as f:
	f.write( '\\documentclass[11pt]{article}\n\n\n' )
	f.write( '\\begin{document}\n\n\n' )
	# write frontmatter:
	for s in frontmatter:
		f.write(  f'{s}\n'  )
	f.write('\n\n\n')
	# write table start:
	f.write(  '\\begin{table}[ht]\n'  )
	f.write(  '\\caption{My caption.}\n'  )
	f.write(  '\\begin{tabular}{ | c | c | c | c | c | c |}\n'  )
	f.write(  '\\hline\n')
	f.write(  'Dataset & Group A & Group B & Mean (B-A) & t value & Significant?\\\\\n')
	f.write(  '\\hline\n')
	# write table rows:
	for i,(y0,y1,md,t,fp) in enumerate( results ):
		s0  = ', '.join( ['%.2f'%yy for yy in y0] )
		s1  = ', '.join( ['%.2f'%yy for yy in y1] )
		mds = '%.3f' %md
		ts  = '%.3f' %t
		fps = 'Yes' if fp else ''
		f.write(  f'{i+1} & \\tiny {s0} & \\tiny {s1} & \\tiny {mds} & \\tiny {ts} & {fps}\\\\\n')
	# write table end:
	f.write(  '\\hline\n')
	f.write(  '\\end{tabular}\n'  )
	f.write(  '\\label{table:sim0d}')
	f.write(  '\\end{table}\n\n\n'  )
	f.write('\\end{document}')








