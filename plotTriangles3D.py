from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import PathPatch
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.pylab as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_toolkits.mplot3d.art3d as art3d
import math
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

filenames = []

filenames.append("output/Bispectrum/Triangle_plots/SN_sparse/LISW_SN_triangle_l100_gaps3.dat")
#filenames.append("output/Bispectrum/Triangle_plots/SN_sparse/LISW_SN_triangle_l200_gaps3.dat")
#filenames.append("output/Bispectrum/Triangle_plots/SN_sparse/LISW_SN_triangle_l300_gaps3.dat")

#zs = [100,200,300]

#for i in [0,1,2]:

f = file(filenames[0])
data_vals = np.loadtxt(f)

    #inverting triangle
#data_vals = np.array([[1,1,20,1,1],[0,1,10,1,0],[0,0,5,0,0]])

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

#ax.annotate('$B_{max} = %s$' % max_val, xy = (0.8,0.52), xytext=(0.8, 0.52), fontsize = 16)

data_reduced = np.zeros((rows,columns))
data = []
for i in range(0,rows):
    for j in range(0,columns):
        if data_new[i][j] == 0:
            data_reduced[i][j] = -1
        else:
            data_reduced[i][j] = math.log10(data_new[i][j]/max_val)
            data.append(data_reduced[i][j])

xs = np.arange(0,columns,1)
ys = np.arange(0,rows,1)
print xs
print ys
xx,yy = np.meshgrid(np.linspace(0,1,columns+1), np.linspace(0,1,rows+1))
X = xx
Y = yy
Z =  10*np.ones(X.shape)
#print Z
for i in range(0,rows):
    for j in range(0,columns):
        if (data_new[i][j] == 0):
            Z[i][j] = np.nan
#print data_reduced
#print Z
#print X
#print Y

norm = colors.LogNorm()
ax.plot_surface(X,Y,Z, rstride=1,cstride=1,facecolors=cm.jet(data_reduced),norm=norm)


plt.show()










#print data_reduced[0]
#print data_new[1]
#print data_reduced[imax][jmax]
#extent = (0,1,0.5,1)
#plt.imshow(data_reduced,extent=extent,interpolation='none')
#divider = make_axes_locatable(ax)

#cax = divider.append_axes("right", size="5%", pad=0.05)
#norm = colors.LogNorm()
#plt.pcolor(data_reduced,norm=norm)
#plt.colorbar(cax=cax)
#plt.title('Small difference between \n $B_{l_1l_2l_3}$ and $B_{max}$', y = 1.03, fontsize=16)
#plt.xlabel('\n \n Large difference between \n $B_{l_1l_2l_3}$ and $B_{max}$', fontsize=16)
#plt.ylabel(r'$\log_{10}\left(\frac{B_{l_1l_2l_3}}{B_{max}}\right)$', fontsize=16)

#plt.show()

