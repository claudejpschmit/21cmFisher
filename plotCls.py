import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size':14})

fig, ax = plt.subplots()

xs1 = []
ys1 = []
xs2 = []
ys2 = []
xs3 = []
ys3 = []
xs4 = []
ys4 = []
xs5 = []
ys5 = []

filename1 = "Cls_140_logscale.dat"
filename2 = "Cl_FG_extragalactic_ff.dat"
filename3 = "Cl_FG_extragalactic_ps.dat"
filename4 = "Cl_FG_gal_ff.dat"
filename5 = "Cl_FG_gal_synch.dat"
f1 = file(filename1)
f2 = file(filename2)
f3 = file(filename3)
f4 = file(filename4)
f5 = file(filename5)
data1 = np.loadtxt(f1)
data2 = np.loadtxt(f2)
data3 = np.loadtxt(f3)
data4 = np.loadtxt(f4)
data5 = np.loadtxt(f5)

ax.set_xscale('log')
ax.set_yscale('log')

plt.xlim([10,10000])
plt.ylim([10E-3,10E8])

print data1[0][0]
print data1[0][1]
for i in range(0, len(data1)):
    xs1.append(data1[i][0])
    ys1.append(data1[i][1])
for i in range(0, len(data2)):
    xs2.append(data2[i][0])
    ys2.append(data2[i][1])
for i in range(0, len(data3)):
    xs3.append(data3[i][0])
    ys3.append(data3[i][1])
for i in range(0, len(data4)):
    xs4.append(data4[i][0])
    ys4.append(data4[i][1])
for i in range(0, len(data5)):
    xs5.append(data5[i][0])
    ys5.append(data5[i][1])

plt.plot(xs1,ys1, label = '21cm Power Spectrum', linewidth = 2)
plt.plot(xs2,ys2, label = 'Extragalactic free-free', linewidth = 2, linestyle = ':')
plt.plot(xs3,ys3, label = 'Extragalactic Point Sources', linewidth = 2, linestyle = '-.')
plt.plot(xs4,ys4, label = 'Galactic free-free', linewidth = 2, linestyle = '--')
plt.plot(xs5,ys5, label = 'Galactic synchrotron', linewidth = 2, linestyle = '-')
plt.legend(loc = 4)
plt.xlabel(r'$l$', fontsize=18)
plt.ylabel(r'$l(l+1)C_l/(2\pi)$ (mK$^2$)', fontsize=18)
plt.title(r'Angular power and foregrounds at $140 MHz$ $(z = 9.2)$')
plt.show()
