import numpy as np
import matplotlib.pyplot as plt
from glob import glob

fig, ax = plt.subplots()

xs = []
ys = []

zmin = 0
zmax = 1000
pkz_steps = 1000
stepsize = (zmax - zmin)/float(pkz_steps)
print stepsize
for n in range(0,pkz_steps):
    ### read data from files...
    z = zmin + n*stepsize
    filename = "CAMB/test_matterpower_" + str(n+1) + ".dat" 
    f = file(filename)
    data = np.loadtxt(f)
    xs.append(z)
    ys.append(data[232][1])
    print f 
    print z, data[232][0], data[232][1]

#f = file("integrand_Ql.dat")
#data = np.loadtxt(f)
#xs2 = data[:,0]
#ys2 = data[:,1]

ys = ys[::-1]

plt.plot(xs,ys)
#plt.plot(xs2,ys2)
plt.show()

