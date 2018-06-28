import numpy as np
import matplotlib as mpl
mpl.use('agg')
import pylab as pl

# Optimized plot parameters:
pl.rcParams['figure.figsize']  = 15, 10
pl.rcParams['figure.dpi']      = 300
pl.rcParams['image.cmap']      = 'jet'
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 30
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'

pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad']  = 8
pl.rcParams['xtick.minor.pad']  = 8
pl.rcParams['xtick.color']      = 'k'
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'

pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad']  = 8
pl.rcParams['ytick.minor.pad']  = 8
pl.rcParams['ytick.color']      = 'k'
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'

# Loading the data from text file:
data = np.loadtxt('data.txt')
N    = data[:, 0]
t    = data[:, 1]

pl.loglog(N, t, '-o', label = r'Mat$\times$Mat') # change label appropriately
# pl.loglog(N, t[0] * N**2 / N[0]**2, '--', color = 'black' label = r'$\mathcal{O}(N^2)$')
textstr = 'Tested Using:\nArrayFire v3.6.0 (OpenCL, 64-bit Linux, build 2858662)\nIntel(R) Xeon(R) CPU E5-2682 v4 @ 2.50GHz, 3951 MB'
pl.text(0.03 * (np.max(N) - np.min(N)), 5e-6 * (np.max(t) - np.min(t)), textstr, fontsize=20)
pl.loglog(N, t[0] * N**3 / N[0]**3, '--', color = 'black', label = r'$\mathcal{O}(N^3)$')
pl.xlabel(r'$N$')
pl.ylabel('Time')
pl.legend()
pl.savefig('plot.png', bbox_inches = 'tight')
