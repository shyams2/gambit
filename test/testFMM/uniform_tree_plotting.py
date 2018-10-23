import numpy as np
import matplotlib.patches as patches
import pylab as pl

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

class Box(object):
    def __init__(self, x_left_bot, y_left_bot, w, h):
        self.x_left_bot = x_left_bot;
        self.y_left_bot = y_left_bot;
        self.w          = w;
        self.h          = h;

    def contains(self, x_coords, y_coords):
        return((x_coords>self.x_left_bot) * (x_coords<(self.x_left_bot + self.w)) * \
               (y_coords>self.y_left_bot) * (y_coords<(self.y_left_bot + self.h))
              )

    def number(self, x_coords, y_coords):
        return(np.sum(self.contains(x_coords, y_coords)))

    def subdivide(self):
        B1 = Box(self.x_left_bot, self.y_left_bot, self.w / 2, self.h / 2);
        B2 = Box(self.x_left_bot +  self.w / 2, self.y_left_bot, self.w / 2, self.h / 2);
        B3 = Box(self.x_left_bot, self.y_left_bot + self.h / 2, self.w / 2, self.h / 2);
        B4 = Box(self.x_left_bot +  self.w / 2, self.y_left_bot +  self.h / 2, self.w / 2, self.h / 2);
        return [B1, B2, B3, B4]

def plot_graph(list_of_boxes):
    fig = pl.figure()
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    ax.set_aspect('equal')
    for B in list_of_boxes[-1]:
        ax.add_patch(patches.Rectangle((B.x_left_bot, B.y_left_bot), B.w, B.h, fill = False))

N        = 5000
x_coords = np.random.randn(N)
y_coords = np.random.randn(N)

# Mapping onto [0, 1]:
x_coords = (x_coords - x_coords.min()) / (x_coords.max() - x_coords.min())
y_coords = (y_coords - y_coords.min()) / (y_coords.max() - y_coords.min())

B_master = Box(0, 0, 1, 1)

n_crit = 10

list_of_boxes = [[B_master]]
def add_level(list_of_boxes):
    temp = []
    for i in range(len(list_of_boxes[-1])):
        temp = temp + list_of_boxes[-1][i].subdivide()
    list_of_boxes.append(temp)

# Checking if reached leaf:
reached_leaf = False
while(reached_leaf == False):
    for i in range(len(list_of_boxes[-1])):
        reached_leaf = True
        if(list_of_boxes[-1][i].number(x_coords, y_coords) > n_crit):
            add_level(list_of_boxes)
            reached_leaf = False

plot_graph(list_of_boxes)
pl.savefig('plot.png', bbox_inches = 'tight')
