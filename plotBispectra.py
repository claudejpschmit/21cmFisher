import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.fftpack
matplotlib.rcParams.update({'font.size':14})

fig, ax = plt.subplots()
xs1 = []
ys1 = []
xs2 = []
ys2 = []
xs3 = []
ys3 = []

filename2 = "Bispectrum_PNG.dat"
filename1 = "LISW_compare.dat"
filename3 = "Bispectrum_equilateral_new.dat"
f1 = file(filename1)
f2 = file(filename2)
f3 = file(filename3)
data1 = np.loadtxt(f1)
data2 = np.loadtxt(f2)
data3 = np.loadtxt(f3)

ax.set_xscale('log')
ax.set_yscale('log')

plt.xlim([2,5000])
plt.ylim([10E-24,10E-11])

for i in range(0, len(data1)):
    xs1.append(data1[i][0])
    ys1.append(data1[i][1])
for i in range(0, len(data2)):
    xs2.append(data2[i][0])
    ys2.append(data2[i][1]/10.0)
for i in range(0, len(data3)):
    xs3.append(data3[i][0])
    ys3.append(data3[i][1]/2.0)

#running average smoothing.
y2 = []
y2.append(ys3[0])
y2.append(ys3[1])
y2.append(ys3[2])
y2.append((ys3[2]+ys3[3]+ys3[4])/3.0)
y2.append((ys3[3]+ys3[4]+ys3[5])/3.0)
y2.append((ys3[4]+ys3[5]+ys3[6])/3.0)
y2.append((ys3[5]+ys3[6]+ys3[7])/3.0)
for i in np.arange(7, len(ys3)-7):
    av = ys3[i-5] +ys3[i-4] + ys3[i-3] + ys3[i-2] + ys3[i-1] + ys3[i] + ys3[i+1] + ys3[i+2] + ys3[i+3] + ys3[i+4] + ys3[i+5]
    av += ys3[i-6] +ys3[i+6]
    av += ys3[i-7] +ys3[i+7]

    av = av / 15.0
    y2.append(av)
 
y2.append(ys3[139])   
y2.append(ys3[140])
y2.append(ys3[141])
y2.append(ys3[142])
y2.append(ys3[143])
y2.append(ys3[144])
y2.append(ys3[145])
plt.plot(xs1,ys1, label = 'LISW', linewidth = 2)
plt.plot(xs2,ys2, label = 'PNG', linewidth = 2, linestyle = ':')
plt.plot(xs3,y2, label = 'Non-linear Gravity', linewidth = 2, linestyle = '--')
#plt.plot(xs3,ys3, label = 'Non-linear Gravity', linewidth = 2)
plt.legend(loc = 3)
plt.xlabel(r'$l$')
plt.ylabel(r'$|B_{lll}|$ (mK$^3$)')
plt.show()
