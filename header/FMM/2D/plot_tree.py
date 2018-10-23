import numpy as np
import h5py
import matplotlib.patches as patches
import pylab as pl

# Setting plot parameters for nicer plots:
pl.rcParams['figure.figsize']  = 9, 4
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

def plot_graph(cx, cy, rx, ry):
    
    fig = pl.figure()

    x_min = (cx - rx).min()
    y_min = (cy - ry).min()
    x_max = (cx + rx).max()
    y_max = (cy + ry).max()

    # For some reason the plot functionality wants
    # things in the range of [0, 1]
    # Mapping values to [0, 1]:
    cx = (cx - x_min) / (x_max - x_min)
    rx = rx / (x_max - x_min)

    cy = (cy - y_min) / (y_max - y_min)
    ry = ry / (y_max - y_min)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    ax.set_aspect('equal')
    
    # Asserting that they are all of the same length:
    try:
        assert(len(cx) == len(cy) == len(rx) == len(ry))
    except:
        raise AssertionError('Elements in array are not the same!!')

    for i in range(len(cx)):
        ax.add_patch(patches.Rectangle((cx[i] - rx[i], cy[i] - ry[i]), 2 * rx[i], 2 * ry[i], fill = False))

h5f = h5py.File('tree_data.h5', 'r')
cx  = h5f['cx'][:]
cy  = h5f['cy'][:]
rx  = h5f['rx'][:]
ry  = h5f['ry'][:]
h5f.close()

plot_graph(cx, cy, rx, ry)
pl.savefig('plot.png', bbox_inches = 'tight')
