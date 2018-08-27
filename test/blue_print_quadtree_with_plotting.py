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
        return(len(self.contains(x_coords, y_coords)))

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

    pl.savefig('plot.png', bbox_inches = 'tight')

x_coords = np.random.rand(100)
y_coords = np.random.rand(100)
B_master = Box(0, 0, 1, 1)

list_of_boxes  = [B_master]
list_of_boxes += B_master.subdivide()
list_of_boxes += list_of_boxes[1].subdivide()
list_of_boxes += list_of_boxes[4].subdivide()
list_of_boxes += list_of_boxes[-1].subdivide()
list_of_boxes += list_of_boxes[-1].subdivide()
list_of_boxes += list_of_boxes[10].subdivide()
list_of_boxes += list_of_boxes[-1].subdivide()
plot_graph(list_of_boxes)
