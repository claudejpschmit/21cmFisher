import matplotlib.pyplot as plt
import numpy as np
import sys

fig, ax = plt.subplots()
filename = sys.argv[1]
filename2 = sys.argv[2]
fisher_mode = sys.argv[3]
data_array = np.loadtxt(filename)
parameters = np.loadtxt(filename2, dtype = 'string')
formatted_params = {'ombh2':r'$\Omega_b h^2$',\
        'omch2': r'$\Omega_{CDM} h^2$',\
        'omega_lambda': r'$\Omega_\Lambda$',\
        'alpha': r'$\alpha$',\
        'beta': r'$\beta$',\
        'gamma': r'$\gamma$',\
        'RLy': r'$R_{Ly}$',\
        'n_s': r'$n_s$',\
        'extragal_ps_A': r'$A^{EG_{ps}}$',\
        'extragal_ps_beta': r'$\beta^{EG_{ps}}$',\
        'extragal_ps_alpha': r'$\alpha^{EG_{ps}}$',\
        'extragal_ps_xi': r'$\xi^{EG_{ps}}$',\
        'extragal_ff_A': r'$A^{EG_{ff}}$',\
        'extragal_ff_beta': r'$\beta^{EG_{ff}}$',\
        'extragal_ff_alpha': r'$\alpha^{EG_{ff}}$',\
        'extragal_ff_xi': r'$\xi^{EG_{ff}}$',\
        'gal_synch_A': r'$A^{G_{synch}}$',\
        'gal_synch_beta': r'$\beta^{G_{synch}}$',\
        'gal_synch_alpha': r'$\alpha^{G_{synch}}$',\
        'gal_synch_xi': r'$\xi^{G_{synch}}$',\
        'gal_ff_A': r'$A^{G_{ff}}$',\
        'gal_ff_beta': r'$\beta^{G_{ff}}$',\
        'gal_ff_alpha': r'$\alpha^{G_{ff}}$',\
        'gal_ff_xi': r'$\xi^{G_{ff}}$'}


plt.imshow(data_array, cmap = 'bone', interpolation = 'nearest')
plt.colorbar()
x = np.arange(0,len(parameters))

labels = [formatted_params[z] for z in parameters]
ax.xaxis.tick_top()
plt.xticks(x, labels, rotation = 90, fontsize = 12)
plt.yticks(x, labels, fontsize = 12)
if fisher_mode == 'inverse':
    plt.title("Inverse Fisher Matrix", y=1.1)
elif fisher_mode == 'fisher':
    plt.title("Fisher Matrix", y = 1.1)
plt.savefig("Matrix_representation.png")
plt.show()
