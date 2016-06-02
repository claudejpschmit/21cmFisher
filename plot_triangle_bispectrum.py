import numpy as np
import copy
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots()
plt.xlabel(r'$l_2/l_1$',fontsize=16)
plt.ylabel(r'$l_3/l_1$', fontsize=16)
ax.annotate(r'folded', xy = (0.5,0.5), xytext=(0.7,0.6), arrowprops=dict(arrowstyle='->',facecolor='black'))
ax.annotate(r'squeezed', xy = (0.0,1.0), xytext=(0.05,0.75), arrowprops=dict(arrowstyle='->',facecolor='black'))
ax.annotate(r'equilateral', xy = (1.0,1.0), xytext=(0.9,0.75), arrowprops=dict(arrowstyle='->',facecolor='black'))

ax.annotate('', xy = (0.01,0.98), xytext=(0.48, 0.51), arrowprops=dict(arrowstyle='<->',facecolor='black'))
ax.annotate('', xy = (0.99,0.98), xytext=(0.52, 0.51), arrowprops=dict(arrowstyle='<->',facecolor='black'))
ax.annotate('elongated', xy = (0.22,0.72),xytext=(0.22, 0.72), rotation=-45)
ax.annotate('isoceles', xy = (0.22,0.72),xytext=(0.72, 0.71), rotation=45)
ax.annotate('$l_1 = 100$', xy = (0.8,0.54),xytext=(0.8, 0.54))
xs1 = []
ys1 = []

filename1 = "output/Bispectrum/Triangle_plots/LISW_triangle_xy.dat"
filename2 = "output/Bispectrum/Triangle_plots/LISW_triangle_values.dat"

f1 = file(filename1)
f2 = file(filename2)
data_xy = np.loadtxt(f1)
data_vals = np.loadtxt(f2)

for i in range(0, len(data_xy)):
    xs1.append(data_xy[i][0])
    ys1.append(data_xy[i][1])

#inverting triangle
data_new = []
rows = len(data_vals)
columns = len(data_vals[0])
for i in range(1,rows):
    data_new.append(data_vals[rows-i])

max_val = 0;
for i in range(0,rows-1):
    for j in range(0,columns-1):
        if max_val < data_new[i][j]:
            max_val = data_new[i][j]

ax.annotate('$B_{max} = %s$' % max_val, xy = (0.8,0.52), xytext=(0.8, 0.52))

data_reduced = copy.copy(data_new)
for i in range(0,rows-1):
    for j in range(0,columns-1):
        if data_new[i][j] == 0:
            data_reduced[i][j] = np.nan
        else:
            data_reduced[i][j] = math.log(data_new[i][j]/max_val)
#print data_reduced
middle = 0.5;
extent = (0,1,0.5,1)
plt.imshow(data_new,extent=extent,interpolation='none')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
plt.title('Small difference between \n $B_{l_1,l_2,l_3}$ and $B_{max}$', y = 1.03, fontsize=16)
plt.xlabel('\n \n Large difference between \n $B_{l_1,l_2,l_3}$ and $B_{max}$')
#ax.annotate('Large difference between \n $B_{l_1,l_2,l_3}$ and $B_{max}$', xy = (0.22,0.72),xytext=(0.94, 0.44), fontsize=16)
plt.ylabel(r'$\ln\left(\frac{B_{l_1,l_2,l_3}}{B_{max}}\right)$', fontsize=16)

plt.show()
