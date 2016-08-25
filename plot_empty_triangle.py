import numpy as np
import copy
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from matplotlib import collections as mc
matplotlib.rcParams.update({'font.size':26})

fig, ax = plt.subplots()
plt.xlabel(r'$l_3/l_1$')
plt.ylabel(r'$l_2/l_1$')
ax.annotate(r'$(c)$', xy = (0.462,0.55))
ax.annotate(r'$(a)$', xy = (0.0,1.0), xytext=(0.05,0.95))
ax.annotate(r'$(e)$', xy = (1.0,1.0), xytext=(0.87,0.95))

ax.annotate('', xy = (0.01,0.98), xytext=(0.48, 0.51), arrowprops=dict(arrowstyle='<->',facecolor='black'))
ax.annotate('', xy = (0.99,0.98), xytext=(0.52, 0.51), arrowprops=dict(arrowstyle='<->',facecolor='black'))
ax.annotate('elongated', xy = (0.2,0.72),xytext=(0.175, 0.72), rotation=-55)
ax.annotate('isoceles', xy = (0.1,0.72),xytext=(0.64, 0.71), rotation=54)

lines = [[(0,1),(0.5,0.5)],[(0.5,0.5),(1,1)]]
c = [(0,0,0),(0,0,0)]
lc = mc.LineCollection(lines, colors=c)
ax.add_collection(lc)

extent = (0,1,0.5,1)
divider = make_axes_locatable(ax)
plt.ylim(0.5,1)
fig.set_size_inches(5,5)
plt.show()
