import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size':14})

fig, ax = plt.subplots()

def cl_fg(l, nu):
    epsilon = 0.000001;
    A1 = 57;
    A2 = 0.014;
    A3 = 700;
    A4 = 0.088;
    n1 = 1.1;
    n2 = 1.0;
    n3 = 2.4;
    n4 = 3.0;
    m1 = 2.07;
    m2 = 2.1;
    m3 = 2.8;
    m4 = 2.15;
    lf = 1000.0;
    nuf = 130.0;

    res = A1*(lf/l)**n1 * (nuf/nu)**m1;
    res += A2*(lf/l)**n2 * (nuf/nu)**m2;
    res += A3*(lf/l)**n3 * (nuf/nu)**m3;
    res += A4*(lf/l)**n4 * (nuf/nu)**m4;
    return epsilon*epsilon*res;


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
xs7 = []
ys7 = []
xs8 = []
ys8 = []
xs9 = []
ys9 = []
xs10 = []
ys10 = []

xs11 = []
ys11 = []

#filename1 = "Cl_400_IM.dat"#"Cls_140_logscale.dat"
#filename2 = "Cl_600_IM.dat"#"Cl_FG_extragalactic_ff.dat"
#filename3 = "Cl_800_IM.dat"#"Cl_FG_extragalactic_ps.dat"
#filename4 = "Cl_1000_IM.dat"#"Cl_FG_gal_ff.dat"
#filename5 = "Cl_1200_IM.dat"#"Cl_FG_gal_synch.dat"

filename1 = "plots/data/test_Cls_z08.dat"#"Cls_140_logscale.dat"
filename2 = "plots/data/test_Cls_z1.dat"#"Cl_FG_extragalactic_ff.dat"
filename3 = "plots/data/test_Cls_z15.dat"#"Cl_FG_extragalactic_ps.dat"
filename4 = "plots/data/test_Cls_z2.dat"#"Cl_FG_gal_ff.dat"
filename5 = "plots/data/test_Cls_z25.dat"#"Cl_FG_gal_synch.dat"

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
filename6 = "plots/data/test_Cls_z08_delta.dat"#"Cls_140_logscale.dat"
filename7 = "plots/data/test_Cls_z1_delta.dat"#"Cl_FG_extragalactic_ff.dat"
filename8 = "plots/data/test_Cls_z15_delta.dat"#"Cl_FG_extragalactic_ps.dat"
filename9 = "plots/data/test_Cls_z2_delta.dat"#"Cl_FG_gal_ff.dat"
filename10 = "plots/data/test_Cls_z25_delta.dat"#"Cl_FG_gal_synch.dat"
filename11 = "plots/data/test_Cl_Noise.dat"#"Cl_FG_gal_synch.dat"

f6 = file(filename6)
f7 = file(filename7)
f8 = file(filename8)
f9 = file(filename9)
f10 = file(filename10)
f11 = file(filename11)
data6 = np.loadtxt(f6)
data7 = np.loadtxt(f7)
data8 = np.loadtxt(f8)
data9 = np.loadtxt(f9)
data10 = np.loadtxt(f10)
data11=np.loadtxt(f11)
ax.set_xscale('log')
ax.set_yscale('log')

#plt.xlim([10,10000])
#plt.ylim([10E-3,10E8])

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
for i in range(0, len(data6)):
    xs6.append(data6[i][0])
    ys6.append(data6[i][1])
for i in range(0, len(data7)):
    xs7.append(data7[i][0])
    ys7.append(data7[i][1])
for i in range(0, len(data8)):
    xs8.append(data8[i][0])
    ys8.append(data8[i][1])
for i in range(0, len(data9)):
    xs9.append(data9[i][0])
    ys9.append(data9[i][1])
for i in range(0, len(data10)):
    xs10.append(data10[i][0])
    ys10.append(data10[i][1])
for i in range(0, len(data11)):
    xs11.append(data11[i][0])
    ys11.append(data11[i][1])

plt.plot(xs1,ys1, label = r'$\nu$ = 400 MHz', linewidth = 2)
plt.plot(xs2,ys2, label = r'$\nu$ = 600 MHz', linewidth = 2)
plt.plot(xs3,ys3, label = r'$\nu$ = 800 MHz', linewidth = 2)
plt.plot(xs4,ys4, label = r'$\nu$ = 1000 MHz', linewidth = 2)
plt.plot(xs5,ys5, label = r'$\nu$ = 1200 MHz', linewidth = 2)
#for i in range(0, len(data1)):
#    xs1.append(data1[i][0])
#    ys1.append(abs(data1[i][1] - data6[i][1]))
#for i in range(0, len(data2)):
#    xs2.append(data2[i][0])
#    ys2.append(abs(data2[i][1] - data7[i][1]))
#for i in range(0, len(data3)):
#    xs3.append(data3[i][0])
#    ys3.append(abs(data3[i][1] - data8[i][1]))
#for i in range(0, len(data4)):
#    xs4.append(data4[i][0])
#    ys4.append(abs(data4[i][1] - data9[i][1]))
#for i in range(0, len(data5)):
#    xs5.append(data5[i][0])
#    ys5.append(abs(data5[i][1] - data10[i][1]))

#print ys5
#print data10
#plt.plot(xs1,ys1, label = r'$\nu$ = 790 MHz, z = 0.8', linewidth = 2)
#plt.plot(xs2,ys2, label = r'$\nu$ = 710 MHz, z = 1', linewidth = 2)
#plt.plot(xs3,ys3, label = r'$\nu$ = 570 MHz, z = 1.5', linewidth = 2)
#plt.plot(xs4,ys4, label = r'$\nu$ = 475 MHz, z = 2', linewidth = 2)
#plt.plot(xs5,ys5, label = r'$\nu$ = 405 MHz, z = 2.5', linewidth = 2)
#plt.plot(xs6,ys6, label = r'$\nu$ = 790 MHz, z = 0.8', linewidth = 2)
#plt.plot(xs7,ys7, label = r'$\nu$ = 710 MHz, z = 1', linewidth = 2)
#plt.plot(xs8,ys8, label = r'$\nu$ = 570 MHz, z = 1.5', linewidth = 2)
#plt.plot(xs9,ys9, label = r'$\nu$ = 475 MHz, z = 2', linewidth = 2)
#plt.plot(xs10,ys10, label = r'$\nu$ = 405 MHz, z = 2.5', linewidth = 2)
plt.plot(xs11,ys11, label = r'Noise, z = 1', linewidth = 2)

ls = []
ys = []
for l in range(1, 10000):
    ls.append(l)
    ys.append(l*(l+1)*cl_fg(l, 710)*1000000/(2*3.1415))


plt.plot(ls,ys, label = r'Foreground residuals, z = 1', linewidth = 2)


plt.legend(loc = 4)
plt.ylim([1e-3,1e5])
plt.xlabel(r'$\ell$', fontsize=18)
plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$ ($\mu$K$^2$)', fontsize=18)
#plt.title(r'Angular power and foregrounds at $140 MHz$ $(z = 9.2)$')
plt.show()
