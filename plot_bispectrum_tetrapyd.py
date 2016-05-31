import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


fig =  plt.figure()
ax = fig.gca(projection='3d')
ys1 = []

filename1 = "LISW_tetrapyd.dat"

f1 = file(filename1)
data = np.loadtxt(f1)

a = [1,2,3]
b = [10,20,30]
a,b = np.meshgrid(a,b)
print a
print b

X = []
Y = []
Z = []
val = []
for i in np.arange(0,len(data)):
    X.append(data[i][0])
    Y.append(data[i][1])
    Z.append(data[i][2])
    val.append(data[i][3])
print Y
X,Y = np.meshgrid(X,Y)
print Y
R = X+Y
surf = ax.plot_surface(
    X, Y, Z, rstride=1, cstride=1,
    facecolors=cm.jet(R),
    linewidth=0, antialiased=False, shade=False)

plt.show()

