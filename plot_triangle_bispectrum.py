import numpy as np
import copy
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.rcParams.update({'font.size':24})
import sys

fig, ax = plt.subplots()
plt.xlabel(r'$\ell_3/\ell_1$')#,fontsize=16)
plt.ylabel(r'$\ell_2/\ell_1$')#, fontsize=16)
ax.annotate(r'folded', xy = (0.5,0.5), xytext=(0.7,0.6), arrowprops=dict(arrowstyle='->',facecolor='black'))
ax.annotate(r'squeezed', xy = (0.0,1.0), xytext=(0.05,0.75), arrowprops=dict(arrowstyle='->',facecolor='black'))
ax.annotate(r'equilateral', xy = (1.0,1.0), xytext=(0.8,0.75), arrowprops=dict(arrowstyle='->',facecolor='black'))

ax.annotate('', xy = (0.01,0.98), xytext=(0.48, 0.51), arrowprops=dict(arrowstyle='<->',facecolor='black'))
ax.annotate('', xy = (0.99,0.98), xytext=(0.52, 0.51), arrowprops=dict(arrowstyle='<->',facecolor='black'))
ax.annotate('elongated', xy = (0.22,0.72),xytext=(0.22, 0.72), rotation=-45)
ax.annotate('isoceles', xy = (0.1,0.72),xytext=(0.69, 0.71), rotation=45)
#ax.annotate('$l_1 = 40$', xy = (0.8,0.54),xytext=(0.8, 0.54))#, fontsize = 16)
#plt.title('Bispectrum due to Non-Linear Gravity')
#plt.title('Bispectrum due to LISW')
#plt.title('Bispectrum due to PNG')

#filename = "output/Bispectrum/Triangle_plots/Bispectrum_PNG_triangle_l40.dat"
#filename = "output/Bispectrum/Triangle_plots/Bispectrum_NLG_triangle_l40.dat"
#filename = "output/Bispectrum/Triangle_plots/LISW_triangle_values_l50.dat"
#filename = "output/Bispectrum/Triangle_plots/SN_sparse/LISW_SN_triangle_l180_gaps0.dat"
#filename = "output/Bispectrum/Triangle_plots/SN_new/LISW_test_triangle_l"
filename = "output/Bispectrum/Triangle_plots/NLG_2/NLG_signal_triangle_l"
number = sys.argv[1]
filename = filename + str(number) + ".dat"
print filename


f = file(filename)
data_vals = np.loadtxt(f)

#inverting triangle
data_new = []
rows = len(data_vals)
columns = len(data_vals[0])
for i in range(1,rows+1):
    data_new.append(data_vals[rows-i])

max_val = 0
imax = 0
jmax = 0
for i in range(0,rows):
    for j in range(0,columns):
        if max_val < data_new[i][j]:
            max_val = data_new[i][j]
            imax = i
            jmax = j

ax.annotate('$B_{max} = %s$' % max_val, xy = (0.7,0.52), xytext=(0.7, 0.52))#, fontsize = 16)
ax.annotate('$\ell_{max} = %s$' % number, xy = (0.2,0.52), xytext=(0.2, 0.52))#, fontsize = 16)


data_reduced = np.zeros((rows,columns))
data = []
for i in range(0,rows):
    for j in range(0,columns):
        if data_new[i][j] == 0:
            data_reduced[i][j] = np.nan
        else:
            data_reduced[i][j] = math.log10(data_new[i][j]/max_val)
            data.append(data_reduced[i][j])
#print imax, jmax, data_new[imax][jmax], max_val

#print data_reduced[0]
#print data_new[1]
#print data_reduced[imax][jmax]
extent = (0,1,0.5,1)
plt.imshow(data_reduced,extent=extent,interpolation='none')
divider = make_axes_locatable(ax)

cax = divider.append_axes("right", size="5%", pad=0.05)
norm = colors.LogNorm()
#plt.pcolor(data_reduced,norm=norm)
plt.colorbar(cax=cax)
#plt.title('Small difference between \n $B_{l_1l_2l_3}$ and $B_{max}$', y = 1.03)#, fontsize=16)
#plt.xlabel('\n \n Large difference between \n $B_{l_1l_2l_3}$ and $B_{max}$')#, fontsize=16)
plt.ylabel(r'$\log_{10}\left(\frac{B_{\ell_1\ell_2\ell_3}}{B_{max}}\right)$' ,fontsize=30)
fig.set_size_inches(5,5)
plt.show()
