import numpy as np
import copy
import math
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
xs1 = []
ys1 = []

filename1 = "LISW_triangle_xy.dat"
filename2 = "LISW_triangle_values.dat"

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

data_reduced = copy.copy(data_new)
for i in range(0,rows-1):
    for j in range(0,columns-1):
        if data_new[i][j] == 0:
            data_reduced[i][j] = 0
        else:
            data_reduced[i][j] = math.log(max_val/data_new[i][j])
print data_reduced
middle = 0.5;
extent = (0,1,0.5,1)
plt.imshow(data_new,extent=extent,interpolation='none')
plt.colorbar()

plt.show()
