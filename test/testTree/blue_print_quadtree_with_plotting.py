import numpy as np
import matplotlib.patches as patches
import pylab as pl
pl.style.use('prettyplot')

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
    for B in list_of_boxes:
        ax.add_patch(patches.Rectangle((B.x_left_bot, B.y_left_bot), B.w, B.h, fill = False))

N        = 4000
x_coords = np.random.rand(N)
y_coords = np.random.rand(N)
B_master = Box(0, 0, 1, 1)

n_crit = 4

list_of_boxes = [B_master]
i             = 0 # counter
while(1):
    try:
        if(list_of_boxes[i].number(x_coords, y_coords) > n_crit):
            list_of_boxes += list_of_boxes[i].subdivide()
        i += 1
    except:
        break

plot_graph(list_of_boxes)
pl.plot(x_coords, y_coords, 'ro')
pl.savefig('plot.png', bbox_inches = 'tight')
