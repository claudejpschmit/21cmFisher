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


