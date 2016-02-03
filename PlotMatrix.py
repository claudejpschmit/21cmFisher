import matplotlib.pyplot as plt
import numpy as np
import sys


filename = sys.argv[1]
data_array = np.loadtxt(filename)
plt.imshow(data_array, cmap = 'bone')
plt.colorbar()
plt.show()
