import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size':14})

ax = plt.subplot(221)

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
xs6 = []
ys6 = []


filename1 = "test_Cls_z08.dat"#"Cls_140_logscale.dat"
filename2 = "test_Cls_z1.dat"#"Cl_FG_extragalactic_ff.dat"
filename3 = "test_Cls_z15.dat"#"Cl_FG_extragalactic_ps.dat"
filename4 = "test_Cls_z2.dat"#"Cl_FG_gal_ff.dat"
filename5 = "test_Cls_z25.dat"#"Cl_FG_gal_synch.dat"
filename6 = "test_Cl_Noise.dat"#"Cl_FG_gal_synch.dat"

f1 = file(filename1)
f2 = file(filename2)
f3 = file(filename3)
f4 = file(filename4)
f5 = file(filename5)
f6 = file(filename6)
data1 = np.loadtxt(f1)
data2 = np.loadtxt(f2)
data3 = np.loadtxt(f3)
data4 = np.loadtxt(f4)
data5 = np.loadtxt(f5)
data6 = np.loadtxt(f6)

ax.set_xscale('log')
ax.set_yscale('log')

#plt.xlim([10,10000])
#plt.ylim([10E-3,10E8])

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
for i in range(0, len(data6)):
    xs6.append(data6[i][0])
    ys6.append(data6[i][1])

plt.plot(xs1,ys1, label = r'$\nu$ = 800 MHz', linewidth = 2)
plt.plot(xs2,ys2, label = r'$\nu$ = 700 MHz', linewidth = 2)
plt.plot(xs3,ys3, label = r'$\nu$ = 550 MHz', linewidth = 2)
plt.plot(xs4,ys4, label = r'$\nu$ = 470 MHz', linewidth = 2)
plt.plot(xs5,ys5, label = r'$\nu$ = 400 MHz', linewidth = 2)
plt.plot(xs6,ys6, label = r'Thermal noise at $\nu$ = 700 MHz', linewidth = 2)
plt.legend(loc = 4)
plt.xlabel(r'$l$', fontsize=18)
plt.ylabel(r'$l(l+1)C_l/(2\pi)$ (mK$^2$)', fontsize=18)
plt.title(r'Angular power spectra between $z = 0.8$ and $z = 2.5$.')

ax = plt.subplot(222)
filename1 = "test_LISW_bispectrum.dat"#"Cl_FG_gal_ff.dat"
filename2 = "test_NLG_bispectrum.dat"#"Cl_FG_gal_synch.dat"
f1 = file(filename1)
f2 = file(filename2)
data1 = np.loadtxt(f1)
data2 = np.loadtxt(f2)

ax.set_xscale('log')
ax.set_yscale('log')

#plt.xlim([10,10000])
#plt.ylim([10E-3,10E8])

xs1 = []
xs2 = []
ys1 = []
ys2 = []

for i in range(0, len(data1)):
    if data1[i][0] != data1[i-1][0]:
        xs1.append(data1[i][0])
        ys1.append(data1[i][1])
for i in range(0, len(data2)):
    if data2[i][0] != data2[i-1][0]:
        xs2.append(data2[i][0])
        ys2.append(data2[i][1])

plt.plot(xs1,ys1, label = r'LISW', linewidth = 2)
plt.plot(xs2,ys2, label = r'NLG', linewidth = 2)
plt.legend(loc = 3)
plt.xlabel(r'$l$', fontsize=18)
plt.ylabel(r'$|B_{lll}|$ (mK$^3$)', fontsize=18)
plt.title(r'Bispectum contributions for equilateral triangles at $z=1$.')


ax = plt.subplot(223)
filename1 = "SN_20_10000_delta100.dat"#"Cl_FG_gal_ff.dat"
f1 = file(filename1)
data1 = np.loadtxt(f1)

#plt.xlim([10,10000])
#plt.ylim([10E-3,10E8])

xs1 = []
ys1 = []

for i in range(0, len(data1)):
    xs1.append(data1[i][0])
    ys1.append(data1[i][1])

plt.plot(xs1,ys1,linewidth = 2)
plt.legend(loc = 4)
plt.xlabel(r'$l$', fontsize=18)
plt.ylabel(r'$S/N$', fontsize=18)
plt.title(r'Signal to noise ratio for LISW Bispectrum at $z = 1$.')

ax = plt.subplot(224)
filename1 = "test_Qls_ref.dat"#"Cl_FG_gal_ff.dat"
f1 = file(filename1)
data1 = np.loadtxt(f1)
ax.set_xscale('log')
ax.set_yscale('log')

#plt.xlim([10,10000])
#plt.ylim([10E-3,10E8])

xs1 = []
ys1 = []

for i in range(0, len(data1)):
    l = data1[i][0] 
    xs1.append(data1[i][0])
    ys1.append(l*(l+1)*data1[i][1])

plt.plot(xs1,ys1,linewidth = 2)
plt.legend(loc = 4)
plt.xlabel(r'$l$', fontsize=18)
plt.ylabel(r'$l (l+1) Q_l / (2\pi)$ (mK)', fontsize=18)
plt.title(r'$Q_l$ at $z = 1$.')

plt.show()
