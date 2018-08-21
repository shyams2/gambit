import numpy as np
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

N   = np.array([2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128])
err = np.zeros(N.size)

for i in range(N.size):
    ana = np.sort(np.loadtxt('nodes_analytic/nodes'+str(N[i])))
    num = np.loadtxt('nodes_numerical/nodes'+str(N[i]))

    err[i] = np.linalg.norm(ana - num)

print(err)
# Gives output as [1.57009246e-16 1.75541673e-16 4.59436527e-16 2.96348453e-16
#                  2.22044605e-16 8.21724430e-16 1.04903592e-15 2.17991733e-15
#                  1.74628677e-15 4.69939119e-15 5.65444499e-15 4.53137466e-15]

# Plotting the figure:
pl.loglog(N, err, '-o')
pl.xlabel(r'$N$')
pl.ylabel(r'$||\mathrm{ChebNodes}_{\mathrm{analytic}} - \mathrm{ChebNodes}_{\mathrm{numerical}}||_2$')
pl.savefig('plot.png')
