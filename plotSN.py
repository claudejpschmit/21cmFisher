import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.fftpack
matplotlib.rcParams.update({'font.size':16})

fig, ax = plt.subplots()
xs1 = []
ys1 = []
xs2 = []
ys2 = []
xs3 = []
ys3 = []

#filename1 = "SN_all4.dat"
filename1 = "SN_min-2_max-10000_delta-100_z-1.dat"
#filename2 = "LISW_compare.dat"
#filename3 = "Bispectrum_equilateral_new.dat"
f1 = file(filename1)
#f2 = file(filename2)
#f3 = file(filename3)
data1 = np.loadtxt(f1)
#data2 = np.loadtxt(f2)
#data3 = np.loadtxt(f3)

#ax.set_xscale('log')
#ax.set_yscale('log')

#plt.xlim([2,5000])
plt.ylim([1E-6,1E-3])


for i in range(0, len(data1)):
    xs1.append(data1[i][0])
    ys1.append(data1[i][1])
#for i in range(0, len(data2)):
#    xs2.append(data2[i][0])
#    ys2.append(data2[i][1]/10.0)
#for i in range(0, len(data3)):
#    xs3.append(data3[i][0])
#    ys3.append(data3[i][1]/2.0)
ys_ones = []
for i in xs1:
    ys_ones.append(1)
plt.plot(xs1,ys1, linewidth = 2)
ax.set_yscale('log')
ax.set_xscale('log')
plt.plot(xs1,ys_ones, linewidth = 1.5, linestyle = '--')
#plt.plot(xs2,ys2, label = 'PNG', linewidth = 2, linestyle = ':')
#plt.plot(xs3,y2, label = 'Non-linear Gravity', linewidth = 2, linestyle = '--')
#plt.plot(xs3,ys3, label = 'Non-linear Gravity', linewidth = 2)
plt.legend(loc = 2)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$S/N$')
plt.show()
